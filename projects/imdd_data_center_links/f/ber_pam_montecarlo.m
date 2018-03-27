function [ber_count, ber_gauss, OSNRdB, Rx] = ber_pam_montecarlo(mpam, Tx, Fibers, Rx, sim)
%% Calculate BER of M-PAM IM-DD system through montecarlo simulation
% Inputs:
% - mpam: PAM class
% - Tx: struct with transmitter paramters
% - Fibers: array of Fiber classes containing the fibers used in
% transmittion i.e., Fibers = [SMF DCF] or Fibers = [SMF];
% - Rx: struct with receiver parameters
% - sim: struct with simulation parameters

mpam = mpam.adjust_levels(Tx.Ptx, sim.rexdB);
mpam = mpam.norm_levels; % normalizes levels but preserves extinction ratio

dataTX = randi([0 mpam.M-1], 1, sim.Nsymb); % Random sequence 
dataTX([1 end]) = 0; % set first and last to zero to ensure periodicity
xd = mpam.signal(dataTX); % Modulated PAM signal

%% ============================ Preemphasis ===============================
if isfield(sim, 'preemphasis') && sim.preemphasis
    femph = abs(freq_time(sim.Nsymb*sim.ros.txDSP, mpam.Rs*sim.ros.txDSP));
    femph(femph >= sim.preemphRange) = 0;
    preemphasis_filter = 10.^(polyval([-0.0013 0.5846 0], femph/1e9)/20);  % Coefficients were measured in the lab  

    xd = real(ifft(fft(xd).*ifftshift(preemphasis_filter)));
end

%% ================================ DAC ===================================
% Set DAC time offset in order to remove group delay due to pulse shaping. 
% This way the first sample of xt will be the center of the first pulse. 
% This is only important for plotting.
% Note: Driving signal xd must be normalized by Vpi
Tx.DAC.offset = sim.Mct/mpam.pulse_shape.sps*(length(mpam.pulse_shape.h)-1)/2;
xt = Tx.DAC.convert(xd, sim.f/sim.fs); % signal varying from -1 to 1

%% ============================= Modulator ================================
if isfield(Tx, 'Mod') % modulator was provided
    % 
    Etx = Tx.Mod.modulate(Tx.Laser.cw(sim.N, sim.fs), xt, sim.f/sim.fs);
else % use laser
    Etx = Tx.Laser.modulate(xt, sim.fs);
end

Pt = abs(Etx).^2;

% Ensures that transmitted power is at the right level
Etx = Etx*sqrt(Tx.Ptx/mean(Pt));

fprintf('Launched power: %.2f dBm\n', power_meter(Etx));

% Calculate extinction ratio and intensity levels position
Psamp = Pt(1:sim.Mct:end);
Pl = zeros(1, mpam.M);
for k = 1:mpam.M
    Pl(k) = mean(Psamp(dataTX == k-1));
end
Pl = sort(Pl, 'ascend');
rexdB = 10*log10(min(Pl)/max(Pl));
fprintf('Estimated extinction ratio = %.2f dB\n', rexdB)
Pl = Pl - Pl(1);
Pl = (mpam.M-1)*Pl/Pl(end);
fprintf('Normalized optical levels: %s\n', sprintf('%.3f\t', Pl));

%% ========================= Fiber propagation ============================
Erx = Etx;
attdB = 0;
% link_gain = Amp.Gain*Rx.PD.R;
for k = 1:length(Fibers)
    Fiberk = Fibers(k); 
    attdB = attdB + Fiberk.att(Tx.Laser.wavelength)*Fiberk.L/1e3;
    Erx = Fiberk.linear_propagation(Erx, sim.f, Tx.Laser.wavelength); % propagation through kth fiber in Fibers
end

%% ========================= Preamplifier =================================
OSNRdB = Inf; % only meaningful when there's a pre-amplifier
if isfield(sim, 'preAmp') && sim.preAmp % only included if sim.preAmp is true
    disp('- IMPORTANT: Simulation includes optical amplifier!')
    [Erx, OSNRdBtheory] = Rx.OptAmp.amp(Erx, sim.fs);
  
    % Measure OSNR
    Osa = OSA(0.3); % optical spectrum analyser with course enough resolution to measure most of signal power
    [~, OSNRdBmeasured] = Osa.estimate_osnr(Erx, Tx.Laser.wavelength, sim.f, sim.shouldPlot('OSNR'));
    % Note: second parameter of Osa.estimate_osnr is OSNR in 0.1 nm
    % resolution
    
    fprintf('OSNR = %.2f dB (theory)\nOSNR = %.2f dB (measured)\n', OSNRdBtheory, OSNRdBmeasured)
    OSNRdB = OSNRdBtheory;
end

%% ========================== Receiver ====================================
%% Direct detection
% Direct detection and add thermal noise
% PD.detect(Ein: Input electric field, fs: sampling rate of samples in Ein,
% noise statistics {'gaussian', 'no noise'}, N0: one-sided PSD of thermal
% noise)
[PrxdBm, Prx] = power_meter(Erx);
yt = Rx.PD.detect(Erx, sim.fs, 'gaussian', Rx.N0);

fprintf('> Received power: %.2f dBm\n', PrxdBm);

% Gain control
yt = yt - mean(yt);
yt = yt*sqrt(mean(abs(mpam.a).^2)/(sqrt(2)*mean(abs(yt).^2)));

%% ADC
% ADC performs filtering, quantization, and downsampling
% For an ideal ADC, ADC.ENOB = Inf
% Rx.ADC.offset = -1; % time offset when sampling      
timeRefSignal = xt;
[yk, ~, ytf] = Rx.ADC.convert(yt, sim.f/sim.fs, timeRefSignal);

%% =========================== Equalization ===============================
Rx.eq.trainSeq = dataTX; % training sequence
mpam = mpam.unbias;
[yd, Rx.eq] = equalize(Rx.eq, yk, [], mpam, sim, sim.shouldPlot('Equalizer'));

% Symbols to be discard in BER calculation
ndiscard = [1:Rx.eq.Ndiscard(1)+sim.Ndiscard (sim.Nsymb-Rx.eq.Ndiscard(2)-sim.Ndiscard):sim.Nsymb];
ydfull = yd;
yd(ndiscard) = []; 
dataTX(ndiscard) = [];

%% =========================== Detection ==============================
% dataRX = mpam.demod(yd);
[dataRX, mpam] = mpam.demod_sweeping_thresholds(yd, dataTX);
% Note: when performing demodulaton sweepingn thresholds, the decision
% thresholds will be picked to minimize BER

%% Counted BER
[~, ber_count] = biterr(dataRX, dataTX);

%% ====================== Gaussian approximation ===========================
% Note 1: Gaussian approximation is used despite the fact that the dominant
% noise in pre-amplified systems is signal-spontaneous beat noise, which is
% not Gaussian distributed.
% Note 2: This calculation includes noise enhancement due to equalization.

% Noise bandwidth
Hrxeq = abs(Rx.ADC.Filter.H(sim.f/sim.fs).*Rx.eq.Hff(sim.f/(Rx.eq.ros*mpam.Rs))).^2;
Hrxeq = Hrxeq/interp1(sim.f, Hrxeq, 0); % normalize to have unit gain at DC
noiseBW = 0.5*trapz(sim.f, Hrxeq);
Rx.noiseBW = noiseBW;

if isfield(sim, 'preAmp') && sim.preAmp % amplified system: signal-spontaneous beat noise dominant
    BWopt = sim.fs; % optical filter (OF) noise bandwidth = sampling rate, since OF is not included 
    % Not divided by 2 because optical filter is a bandpass filter

    % Noise std for intensity level Plevel
    Npol = 2; % number of polarizations. Npol = 1, if polarizer is present, Npol = 2 otherwise.
    
    noiseSTD = @(Plevel) sqrt(noiseBW*Rx.N0 + Rx.PD.varShot(Plevel, noiseBW)... % thermal + shot
                    + Rx.PD.R^2*(Rx.OptAmp.varNoiseDD(Plevel/(Rx.OptAmp.Gain), noiseBW, BWopt, Npol))); % sig-spont + spont-spont    
    % Note 1: Plevel is divided by amplifier gain to obtain power at the amplifier input
    % Note 2: if the amplifier is operating in the constant output power 
    % mode, the appropriate amplifier gain was already set in montecarlo
    % simulation. The BER estimation function will assume the same amplifier gain.

else % unamplified system: thermal-noise dominant
    noiseSTD = Rx.PD.stdNoise(Rx.ADC.Filter.H(sim.f/sim.fs), Rx.eq.Hff(sim.f/(Rx.eq.ros*mpam.Rs)), Rx.N0, Tx.Laser.RIN, sim);
end

% AWGN approximation
mpamRef = mpam.adjust_levels(Prx, -Inf); % may replace -Inf for rexdB
ber_gauss = mpamRef.berAWGN(noiseSTD);

%% Plots
if sim.shouldPlot('Optical eye diagram')  
    Nstart = sim.Ndiscard*sim.Mct + 1;
    figure(103), clf, hold on, box on
    eyediagram(abs(Etx(Nstart:end)).^2, 2*sim.Mct)
    a = axis;
    plot(a(1:2), [1; 1]*Pl, '-k', 'Linewidth', 2)
    title('Transmitted optical signal eye diagram')
    drawnow
end

if sim.shouldPlot('Received signal eye diagram')
    mpam = mpam.norm_levels();
    Nstart = sim.Ndiscard*sim.Mct + 1;
    figure(104), clf, box on, hold on
    eyediagram(ytf(Nstart:end), 2*sim.Mct)
    title('Received signal eye diagram')
    a = axis;
    h1 = plot(a(1:2), (mpam.a*[1 1]).', '-k');
    h2 = plot(a(1:2), (mpam.b*[1 1]).', '--k');
    h3 = plot((sim.Mct+1)*[1 1], a(3:4), 'k');
    legend([h1(1) h2(1) h3], {'Levels', 'Decision thresholds', 'Sampling point'})
    drawnow
end

if sim.shouldPlot('Signal after equalization')
    mpam = mpam.norm_levels();
    figure(105), clf, box on, hold on
    h1 = plot(ydfull, 'o');
    a = axis;
    h2= plot(a(1:2), (mpam.a*[1 1]).', '-k');
    h3 = plot(a(1:2), (mpam.b*[1 1]).', '--k');
    h4 = plot((Rx.eq.Ndiscard(1)+sim.Ndiscard)*[1 1], a(3:4), ':k');
    plot(((sim.Nsymb-Rx.eq.Ndiscard(2)-sim.Ndiscard))*[1 1], a(3:4), ':k');
    legend([h1 h2(1) h3(1) h4], {'Equalized samples', 'PAM levels',...
        'Decision thresholds', 'BER measurement window'})
    title('Signal after equalization')
%     axis([1 sim.Nsymb -0.2 1.2])
    drawnow
end

if sim.shouldPlot('Heuristic noise pdf')
    figure(106)
    [nn, xx] = hist(yd(dataTX == 2), 50);
    nn = nn/trapz(xx, nn);
    bar(xx, nn)
    drawnow
end

if sim.shouldPlot('Decision errors')
    figure(107)
    stem(dataRX ~= dataTX)
    ylabel('Errors')
    xlabel('Symbol')
    title('Decision errors')
    drawnow
end