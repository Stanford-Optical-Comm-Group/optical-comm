%% Amplified or unamplified M-PAM IM-DD system
% - Equalization is done using a fractionally spaced linear equalizer
% Simulations include modulator, fiber, optical amplifier (optional) characterized 
% only by gain and noise figure, optical bandpass filter, photodiode (pin or apd),
% antialiasing filter, sampling, and linear equalization
clear, clc

addpath ../../components/
addpath ../../measurement/
addpath ../../signal/
addpath ../../utilities/
addpath ../../processing/
addpath f/

%% Transmit power swipe
Tx.PtxdBm = -25:-15; % transmitter power range

%% Simulation parameters
sim.Rb = 112e9;    % bit rate in bits/sec
sim.Nsymb = 2^16; % Number of symbols in montecarlo simulation
sim.ros.txDSP = 1; % oversampling ratio transmitter DSP (must be integer). DAC samping rate is sim.ros.txDSP*mpam.Rs
% For DACless simulation must make Tx.dsp.ros = sim.Mct and DAC.resolution = Inf
sim.ros.rxDSP = 2; % oversampling ratio of receiver DSP. If equalization type is fixed time-domain equalizer, then ros = 1
sim.Mct = 2*5;      % Oversampling ratio to simulate continuous time. Must be integer multiple of sim.ros.txDSP and numerator of sim.ros.rxDSP
sim.BERtarget = 1.8e-4; 
sim.Ndiscard = 512; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation
sim.Modulator = 'DML'; % 'MZM' or 'DML' (DML also works as EML)
sim.rexdB = -10; % extinction ratio
 
%% Simulation control
sim.preAmp = true;
sim.preemphasis = false; % preemphasis to compensate for transmitter bandwidth limitation
sim.preemphRange = 25e9; % preemphasis range
sim.PMD = true; % whether to simulate PMD
sim.stopSimWhenBERReaches0 = true; % stop simulation when counted BER reaches 0

% Control what should be plotted
sim.Plots = containers.Map();
sim.Plots('BER') = 1;
sim.Plots('Received signal eye diagram') = 0;
sim.Plots('Signal after equalization') = 1;
sim.Plots('Equalizer') = 1;
sim.Plots('Electronic predistortion') = 0;
sim.Plots('Channel frequency response') = 0;
sim.Plots('OSNR') = 0;
sim.Plots('Decision errors') = 0;
sim.Plots('Received signal optical spectrum') = 0;
sim.Plots('PAM levels MZM predistortion') = 0;
sim.shouldPlot = @(x) sim.Plots.isKey(x) && sim.Plots(x);

%% Pulse shape
Tx.pulse_shape = select_pulse_shape('rect', sim.ros.txDSP);
% pulse_shape = select_pulse_shape('rrc', sim.ros.txDSP, 0.5, 6);

%% M-PAM
% PAM(M, bit rate, leve spacing : {'equally-spaced', 'optimized'}, pulse
% shape: struct containing properties of pulse shape 
mpam = PAM(4, sim.Rb, 'equally-spaced', Tx.pulse_shape);

%% Time and frequency
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'
[sim.f, sim.t] = freq_time(sim.N, sim.fs);

%% ==========================  Transmitter ================================
%% DAC
% DAC constructor: DAC(Nhold (samples), resolution (bits), Filter (class Filter))
% Nhold : upsampling factor in samples (must be integer)
% resolution : DAC resolution in bits
% Filter : instance of class Filter to realize analog filtering
F = Filter('bessel', 5, 30e9/(sim.fs/2)); % anti-aliasing filter
Tx.DAC = DAC(sim.Mct/sim.ros.txDSP, 5, F);

%% Laser
% Laser constructor: laser(lambda (nm), PdBm (dBm), RIN (dB/Hz), linewidth (Hz), frequency offset (Hz))
% lambda : wavelength (nm)
% PdBm : output power (dBm)
% RIN : relative intensity noise (dB/Hz)
% linewidth : laser linewidth (Hz)
% freqOffset : frequency offset with respect to wavelength (Hz)
Tx.Laser = Laser(1250e-9, 0, -Inf, 0, 0);
% Tx.Laser.Filter = Filter('two-pole', 30e9, sim.fs); % only necessary if laser is modulated
Tx.Laser.Filter = Filter('bessel', 5, 30e9/(sim.fs/2)); % only necessary if laser is modulated

%% Modulator
if strcmpi(sim.Modulator, 'mzm')
    %% Mach Zenhder modulator (optional)
    F = Filter('two-pole', 30e9, sim.fs); % pre-defined BW limitation
    Tx.Mod = MachZehnderMod(F); 
%     Tx.Mod = MachZehnderMod(0.98, 0.1, sim.f); % loss and velocity mismatch 
    Tx.Mod.Vpi = 1; % since driving signal varies from -1 to 1
    Tx.Mod.Vbias = -0.05; % M-PAM signal is already centered around ~0.5
%     Note: if Tx.Mod is not provided apd_ber assumes that the Laser is
%     directly modulated. This can be used to simulate EML modulators as well 
%     as directly modulated lasers 
%     To achieve AWGN performance it is necessary to compensate for MZM
%     distortion
end

%% ============================= Fiber ====================================
% fiber(Length in m, anonymous function for attenuation versus wavelength
% (default: att(lamb) = 0 i.e., no attenuation), anonymous function for 
% dispersion versus wavelength (default: SMF28 with lamb0 = 1310nm, S0 = 0.092 s/(nm^2.km))
SMF = Fiber(0e3); 
DCF = Fiber(0, @(lamb) 0, @(lamb) -0.1*(lamb-1550e-9)*1e3 - 40e-6); 

Fibers = [SMF DCF];

linkAttdB = SMF.att(Tx.Laser.wavelength)*SMF.L/1e3...
    + DCF.att(Tx.Laser.wavelength)*DCF.L/1e3;

%% ========================== Amplifier ===================================
% Constructor: OpticalAmplifier(Operation, param, Fn, Wavelength)
% - Opertation: either 'ConstantOutputPower' or 'ConstantGain'
% - param: GaindB if Operation = 'ConstantGain', or outputPower
% if Operation = 'ConstantOutputPower'
% - Fn:  noise figure in dB
% - Wavelength: operationl wavelength in m
Rx.OptAmp = OpticalAmplifier('ConstantOutputPower', 0, 5, Tx.Laser.wavelength); 
% Rx.OptAmp = OpticalAmplifier('ConstantGain', 20, 5, Tx.Laser.wavelength);

%% ============================ Receiver ==================================
%% Photodiode
% PIN
% pin(R: Responsivity (A/W), Id: Dark current (A), BW: bandwidth (Hz))
% PIN frequency response is modeled as a first-order system
Rx.PD = PIN(1, 10e-9, Inf);

% APD
% constructor: apd(GaindB : gain in dB, ka : impact ionization factor, BW: bandwidth (Hz), R: Responsivity (A/W), Id: Dark current (A))
% Rx.PD = APD(10, 0.18, [24e9 290e9], 0.74, 40e-9);
% Rx.PD = APD(10, 0.18, Inf, 0.74, 40e-9);

%% TIA-AGC
% One-sided thermal noise PSD
Rx.N0 = (30e-12).^2; 

%% Receiver DSP
%% ADC
% ADC constructor: ADC(T (samples), ENOB (bits), Filter (class Filter))
% T : sampling period in samples (does not need to be an integer)
% ENOB : effective number of bits
% Filter : instance of class Filter to realize anti-aliasing filter
fs = sim.ros.rxDSP*mpam.Rs*(mpam.pulse_shape.rolloff + 1)/2; % sampling frequency
F = Filter('butter', 5, 0.5*fs/(sim.fs/2)); % anti-aliasing filter
Rx.ADC = ADC(sim.Mct/sim.ros.rxDSP, 5, F);

%% Equalizer
% Terminology: TD = time domain, SR = symbol-rate, LE = linear equalizer
Rx.eq.ros = sim.ros.rxDSP;
Rx.eq.type = 'Adaptive TD-LE';
Rx.eq.Ntaps = 15; % number of taps
Rx.eq.mu = 5e-3; % adaptation ratio
Rx.eq.Ntrain = 1e4; % Number of symbols used in training (if Inf all symbols are used)
Rx.eq.Ndiscard = [1.2e4 1024]; % symbols to be discard from begining and end of sequence due to adaptation, filter length, etc

%% Generate summary
% check if there are enough symbols to perform simulation
Ndiscarded = sum(Rx.eq.Ndiscard) + 2*sim.Ndiscard;
assert(Ndiscarded < sim.Nsymb, 'There arent enough symbols to perform simulation. Nsymb must be increased or Ndiscard must be reduced')
fprintf('%d (2^%.2f) symbols will be used to calculated the BER\n', sim.Nsymb - Ndiscarded, log2(sim.Nsymb - Ndiscarded));
Nreq_symbs = -log(1 - 0.95)/(log2(mpam.M)*sim.BERtarget);
fprintf('2^%.2f symbols are necessary to measured target BER with confidence level of 95 %%\n', log2(Nreq_symbs)) 

%% Run simulation
%% Calculate BER
Ptx = dBm2Watt(Tx.PtxdBm); % Transmitted power

ber.gauss = zeros(size(Ptx));
ber.count = zeros(size(Ptx));
OSNRdB = zeros(size(Ptx));
for k = 1:length(Ptx)
    Tx.Ptx = Ptx(k);
         
    % Montecarlo simulation
    [ber.count(k), ber.gauss(k), OSNRdB(k), Rx] = ber_pam_montecarlo(mpam, Tx, Fibers, Rx, sim);
    % only equalizer and noiseBW fields are changed in struct Rx
    
    if sim.stopSimWhenBERReaches0 && ber.count(k) == 0
        break;
    end 
end

%% Plots
if sim.shouldPlot('BER') && length(ber.count) > 1
    figure(1), hold on, box on
    if sim.preAmp
        hline(1) = plot(OSNRdB, log10(ber.count), '-o', 'LineWidth', 2, 'DisplayName', 'Counted');
        hline(2) = plot(OSNRdB, log10(ber.gauss), '-', 'Color', get(hline(1), 'Color'), 'LineWidth', 2, 'DisplayName', 'Gaussian Approx.');
        hline(3) = plot(OSNRdB(OSNRdB ~= 0), log10(pam_ber_from_osnr(mpam.M, OSNRdB(OSNRdB ~= 0), Rx.noiseBW)), '--k', 'LineWidth', 2, 'DisplayName', 'Sig-spont & noise enhancement');
        hline(3) = plot(OSNRdB(OSNRdB ~= 0), log10(pam_ber_from_osnr(mpam.M, OSNRdB(OSNRdB ~= 0), mpam.Rs/2)), ':k', 'LineWidth', 2, 'DisplayName', 'Sig-spont limit');
        legend('-DynamicLegend')
        axis([OSNRdB(1) OSNRdB(find(OSNRdB ~= 0, 1, 'last')) -8 0])
        xlabel('OSNR (dB)', 'FontSize', 12)
    else
        PrxdBm = Tx.PtxdBm - linkAttdB;
        hline(1) = plot(PrxdBm, log10(ber.count), '-o', 'LineWidth', 2, 'DisplayName', 'Counted');
        hline(2) = plot(PrxdBm, log10(ber.gauss), '-', 'Color', get(hline(1), 'Color'), 'LineWidth', 2, 'DisplayName', 'Gaussian Approximation');
        legend('-DynamicLegend')
        axis([PrxdBm(1) PrxdBm(end) -8 0])
        xlabel('Received power (dB)', 'FontSize', 12)
    end
end
ylabel('log_{10}(BER)', 'FontSize', 12)
set(gca, 'FontSize', 12)