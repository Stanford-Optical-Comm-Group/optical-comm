classdef Laser
    %% Single-mode laser
    % Modeling includes intensity noise, phase noise, chirp, frequency
    % offset, and frequency response for directly modulated lasers
    properties
        lambda % wavelength (m)
        PdBm   % Average optical power (dBm)
        RIN    % Relative intensity noise (dB/Hz)
        linewidth % Linewidth (Hz)
        freqOffset = 0;  % frequency offset with respect to wavelength lambda (Hz)
        alpha  = 0; % laser chirp
        Filter = Filter() % instance of class Filter corresponding to bandwidth limitation of laser (only used for modulation)
    end
    
    properties(Dependent)
        PW % power in Watts
        wavelength % same as lambda
    end      
    
    methods
        function obj = Laser(lambda, PdBm, RIN, linewidth, freqOffset)
            %% Class constructor
            obj.lambda = lambda;
            obj.PdBm = PdBm;
            
            if exist('RIN', 'var')
                obj.RIN = RIN;
            end
            
            if exist('linewidth', 'var')
                obj.linewidth = linewidth;
            end
            
            if exist('freqshift', 'var')
                obj.freqOffset = freqOffset;
            end         
        end
        
        function LaserTable = summary(self)
            %% Generate table summarizing class values
            disp('-- Laser class parameters summary:')
            rows = {'Wavelength'; sprintf('Maximum frequency offset from %.2f nm', self.lambda*1e9);...
                'Power'; 'Relative intensity noise'; 'Linewidth'};
            Variables = {'lambda'; 'freqOffset'; 'PdBm'; 'RIN'; 'linewidth'};
            Values = [self.lambda*1e9; max(self.freqOffset)/1e9; self.PdBm; self.RIN; self.linewidth/1e3];
            Units = {'nm'; 'GHz'; 'dBm'; 'dB/Hz'; 'kHz'};

            LaserTable = table(Variables, Values, Units, 'RowNames', rows);
        end
        
        %% Main Methos
        function sigma2 = varRIN(self, P, Df)
            %% RIN variance
            % Inputs
            % - P: Laser power (W)
            % - Df: One-sided bandwidth (Hz)
            sigma2 = [];
            if ~isempty(self.RIN)
                sigma2 = 10^(self.RIN/10)*P.^2*Df; 
            end
        end
        
        function sigma2 = varPN(self, fs)
            %% Variance of phase noise
            % Inputs:
            % - fs: sampling frequency
            sigma2 = [];
            if ~isempty(self.linewidth)
                sigma2 = 2*pi*self.linewidth/fs; 
            end
        end
        
        function Pout = addIntensityNoise(self, Pout, fs)
            %% Add intensity noise to signal Pin sampled at fs
            % Inputs: 
            % Pin: optical power
            % fs: sampling rate
            if not(isinf(self.RIN)) && not(isempty(self.RIN))
                Pout = Pout + sqrt(self.varRIN(Pout, fs/2)).*randn(size(Pout));
            end
        end
        
        function [Eout, phase_noise] = addPhaseNoise(self, Eout, fs)
            %% Add phase noise
            % Inputs:
            % - Ein: Electric field
            % - fs: sampling rate
            if self.linewidth ~= 0 && not(isempty(self.linewidth))
                initial_phase    = pi*(2*rand(1)-1); % [-pi, pi]
                dtheta      =       [0 sqrt(self.varPN(fs))*randn(1, length(Eout)-1)]; % i.i.d.Gaussian random variables with zero mean and variance sigma_p^2
                phase_noise = initial_phase + cumsum(dtheta, 2);
                Eout = Eout.*exp(1j*phase_noise); % adds phase noise
            else
                phase_noise = 0;
            end
        end
        
        function Eout = addTransientChirp(self, Ein)
            %% Add transient chirp
            Eout = Ein.*exp(1j*self.alpha/2*log(abs(Ein).^2));
        end
        
        function Eout = cw(self, N, fs)
            %% Generates continous wave waveform
            % Inputs:
            % - N: number of samples
            % - fs: sampling frequency
                       
            Pout = self.PW*ones(1, N);
            
            Pout = self.addIntensityNoise(Pout, fs);
                        
            Eout = sqrt(Pout);
            
            Eout = self.addPhaseNoise(Eout, fs);
                        
            dt = 1/fs;
            t = 0:dt:((N-1)*dt); % time vector
            if length(self.freqOffset) == 1
                if  self.freqOffset ~= 0
                	Eout = Eout.*exp(1j*2*pi*self.freqOffset*t);
                end
            else
                Eout = Eout.*exp(1j*2*pi*self.freqOffset.*t);
            end            
        end    
        
        function Eout = modulate(self, x, fs)
            %% Modulate laser output based on driving signal x
            % This generates an electric signal with average power (PdBm)
            % The laser output is filtered by the frequency response H
            % If after filtering the output signal becomes negative,
            % clipping is done
            % sim specifies sampling rate (fs), total number of points (N),
            % and whether to include intensity noise (RIN) 
            % Inputs:
            %  - x: driving signal
            %  - sim: simulation struct
            Pout = self.PW*x/mean(x); % adjust to laser power
            
            Pout = self.Filter.filter(Pout); % filter
                                  
            Pout = self.addIntensityNoise(Pout, fs); % intensity noise
            
            % Clip
            Pout(Pout < 0) = 0;
            Eout = sqrt(Pout);
                                    
            Eout = self.addTransientChirp(Eout); % transient chirp
            
            Eout = self.addPhaseNoise(Eout, fs); % phase noise
        end          
    end
    
    %% Get and Set methods
    methods 
        function PW = get.PW(self)
            PW = dBm2Watt(self.PdBm);
        end
        
        function self = set.PW(self, Ptx)
            self.PdBm = Watt2dBm(Ptx);
        end
        
        function self = setPower(self, PtxW)
            self.PdBm = Watt2dBm(PtxW);
        end
        
        function lamb = get.wavelength(self)
            lamb = self.lambda;
        end
        
        function self = set.wavelength(self, lamb)
            self.wavelength = lamb;
            self.lambda = lamb;
        end
    end
end
