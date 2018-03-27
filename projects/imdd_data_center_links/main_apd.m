%% APD-based short-reach system 
% BER is computed using Montecarlo simulation, enumeration, and AWGN approximation including noise enhancement
% Equalization can be adaptive or fixed. If fixed equalization, receiver
% filter is an matched filter and symbol-rate sampling is performed.
% Fixed adaptation:
% Photodiode -> Matched filter -> symbol-rate sampling -> equalizer
% - sim.ros.rxDSP = 1; % oversampling ratio of receivre DSP
% - Rx.filtereing = 'matched'
% - Rx.eq.type = 'Fixed TD-SR-LE';
% Adaptive adaptation:
% Photodiode -> Antialiasing filter -> sampling -> adaptive equalizer
% - sim.ros.rxDSP = 2; % or 5/4, etc, oversampling ratio of receivre DSP
% - Rx.filtereing = 'antialiasing'
% - Rx.eq.type = 'Adaptive TD-LE';

clear, clc

addpath ../../components/
addpath ../../measurement/
addpath ../../signal/
addpath ../../utilities/
addpath ../../processing/
addpath f/

%% Simulation parameters
sim.Nsymb = 2^16; % Number of symbols in montecarlo simulation
sim.ros.txDSP = 1; % oversampling ratio of transmitter DSP
sim.ros.rxDSP = 2; % oversampling ratio of receivre DSP
sim.Mct = 8*sim.ros.rxDSP;    % Oversampling ratio to simulate continuous time (must be odd so that sampling is done right, and FIR filters have interger grpdelay)  
sim.L = 4;        % de Bruijin sub-sequence length (ISI symbol length)
sim.BERtarget = 1.8e-4; 
sim.Ndiscard = 256;  % number of 0 symbols to be inserted at the begining and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation
sim.WhiteningFilter = true;
sim.rexdB = -10; % extinction ratio (Pmin/Pmax)

%
sim.terminateWhenBERReaches0 = true; % whether simulation is terminated when counted BER reaches 0

sim.Plots = containers.Map();
sim.Plots('BER') = 1;
sim.Plots('Adaptation MSE') = 0;
sim.Plots('Frequency Response') = 0;
sim.Plots('Equalizer') = 1;
sim.Plots('Optical eye diagram') = 0;
sim.Plots('Received signal eye diagram') = 0;
sim.Plots('Signal after equalization') = 1;
sim.shouldPlot = @(x) sim.Plots.isKey(x) && sim.Plots(x);

%% M-PAM
pulse_shape = select_pulse_shape('rect', sim.ros.txDSP);
mpam = PAM(4, 112e9, 'equally-spaced', pulse_shape);

%% Time and frequency
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'
[sim.f, sim.t] = freq_time(sim.N, sim.fs);

%% Laser
% Laser constructor: laser(lambda (nm), PdBm (dBm), RIN (dB/Hz), linewidth (Hz), frequency offset (Hz))
% lambda : wavelength (nm)
% PdBm : output power (dBm)
% RIN : relative intensity noise (dB/Hz)
% linewidth : laser linewidth (Hz)
% freqOffset : frequency offset with respect to wavelength (Hz)
Tx.Laser = Laser(1250e-9, 0, -Inf, 0, 0);
Tx.Laser.Filter = Filter('two-pole', 30e9, sim.fs); % only necessary if laser is modulated

%% Transmitter
Tx.PtxdBm = -22:-10; % transmitted power sweep
% Tx.PtxdBm = -15; 

%% DAC
% DAC constructor: DAC(Nhold (samples), resolution (bits), Filter (class Filter))
% Nhold : upsampling factor in samples (must be integer)
% resolution : DAC resolution in bits
% Filter : instance of class Filter to realize analog filtering
F = Filter('bessel', 5, 30e9/(sim.fs/2)); % anti-aliasing filter
Tx.DAC = DAC(sim.Mct/sim.ros.txDSP, Inf, F);

%% Mach Zenhder modulator (optional)
% F = Filter('two-pole', 30e9, sim.fs);
% Tx.Mod = MachZehnderMod(F); 
% Tx.Mod.Vpi = 1;
% Tx.Mod.Vbias = -0.05; % M-PAM signal is already centered around ~0.5
% Note: if Tx.Mod is not provided apd_ber assumes that the Laser is
% directly modulated. This can be used to simulate EML modulators as well 
% as directly modulated lasers 
% To achieve AWGN performance it is necessary to compensate for MZM
% distortion

%% Fiber
SMF = Fiber(0e3);

%% Receiver
Rx.N0 = (30e-12).^2; % thermal noise psd
Rx.filtering = 'antialiasing'; % Electric Lowpass Filter prior to sampling and equalization: either 'antialisaing' or 'matched'

%% ADC
% ADC constructor: ADC(T (samples), ENOB (bits), Filter (class Filter))
% T : sampling period in samples (does not need to be an integer)
% ENOB : effective number of bits
% Filter : instance of class Filter to realize anti-aliasing filter
fs = sim.ros.rxDSP*mpam.Rs*(mpam.pulse_shape.rolloff + 1)/2;
F = Filter('butter', 5, 0.5*fs/(sim.fs/2)); % anti-aliasing filter
Rx.ADC = ADC(sim.Mct/sim.ros.rxDSP, Inf, F);

%% Equalization
% Rx.eq.type = 'None';
Rx.eq.type = 'Adaptive TD-LE';
% Rx.eq.type = 'Fixed TD-SR-LE';
Rx.eq.ros = sim.ros.rxDSP;
Rx.eq.Ntaps = 15;
Rx.eq.mu = 5e-3;
Rx.eq.Ntrain = Inf; % Number of symbols used in training (if Inf all symbols are used)
Rx.eq.Ndiscard = [2e4 1024]; % symbols to be discard from begining and end of sequence due to adaptation, filter length, etc

%% APD 
% apd(GaindB, ka, [BW GBP=Inf], R, Id) 
Apd = APD(10, 0.1, [20e9 270e9], 0.74, 10e-9);
% Apd = apd(15, 0.5, Inf, 1, 10e-9);

% BER
sim.OptimizeGain = true;
ber_apd = apd_ber(mpam, Tx, SMF, Apd, Rx, sim);
% 
mpam.level_spacing = 'optimized';
ber_apd_eq = apd_ber(mpam, Tx, SMF, Apd, Rx, sim);
        