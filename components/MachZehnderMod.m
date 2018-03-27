classdef MachZehnderMod 
    %% Mach-Zenhder modulator
    % Filter by "H" and apply nolinear transfer function "sin"
    properties
        Vpi = 1 % driving voltage to cause a pi phase swing in the electric field response
        Vbias = 0 % bias voltage
        Filter % instance of class Filter that models frequency limitation of MZ modulator
        ratio = 0.98 % velocity mismatch ratio
        L = 0.1 % iteractive length
    end
    
    properties (Dependent)
        BW % modulator bandwidth
    end
    
    properties (Constant, Hidden)
        n_r = 2.15;         % refractive index of coplanar-waveguide (CPW) for TM input light
        h = 6.62606957e-34; % Planck
        q = 1.60217657e-19; % electron charge
        c = 299792458;      % speed of light
    end
    
    methods
        function obj = MachZehnderMod(varargin)
            %% Constructor of class MachZehnderMod
            % Input Filter is the an instance of class Filter corresponding
            % to the bandwidth limitation of the MZ modulator
            
            if isa(varargin{1}, 'Filter') % 
                obj.Filter = varargin{1};
            else
                obj.ratio = varargin{1};
                obj.L = varargin{2}; % 
                f = varargin{3}; % frequency vector
                obj.Filter = Filter(obj.mismatch(f, 'false')); 
            end
        end
        
        function BW = get.BW(self)
            %% Modulator bandwidth
            H0 = abs(self.Filter.H(0)).^2;
            [f, ~, exitflag] = fzero(@(f) abs(self.Filter.H(abs(f)*1e9)).^2/H0 - 0.5, 20);
            if exitflag ~= 1
                warning('Bandwidth estimation did not converge')
            end
            BW = abs(f)*1e9;
        end
        
        function Eout = modulate(self, Ein, Vin, f)
            %% Mach-Zehnder modulator
            % - Dual-pol I/Q modulator if Vin is 2 x N and complex
            % - Single-pol I/Q modulator if Vin is 1 x N and complex
            % - Intensity modulator if Vin is 1 x N and real
            % Inputs:
            % - Ein : input electric field
            % - Vin : driving signal normalized by Vpi/2. Input signal must already
            % have proper biasing and voltage swing
            % - tx : contains modulators parameters
            %   - Mod.Hel : electric frequency response
            % Output : Eout is the modulated electric field
            Hmod = ifftshift(self.Filter.H(f));    

            if size(Vin, 1) == 2 && not(isreal(Vin)) 
            %% Dual-pol I/Q
                % Breaks data into in-phase, quadrature in two pols
                Vxi = real(Vin(1, :))/self.Vpi + self.Vbias;
                Vxq = imag(Vin(1, :))/self.Vpi + self.Vbias;
                Vyi = real(Vin(2, :))/self.Vpi + self.Vbias;
                Vyq = imag(Vin(2, :))/self.Vpi + self.Vbias;

                % Filter driving signals
                Vxi = real(ifft(fft(Vxi).*Hmod));
                Vxq = real(ifft(fft(Vxq).*Hmod)); 
                Vyi = real(ifft(fft(Vyi).*Hmod)); 
                Vyq = real(ifft(fft(Vyq).*Hmod));

                % Modulate signal fields (each has unit average power)
                Enorm = 1/sqrt(sqrt(2)); % normalize so that E(|Vout|^2) = 1
                Vout   = Enorm*[sin(pi*Vxi/2) + 1i*sin(pi*Vxq/2);...        % x polarization
                                sin(pi*Vyi/2) + 1i*sin(pi*Vyq/2)];         % y polarization

                Eout =  [Ein/sqrt(2).*Vout(1, :); Ein/sqrt(2).*Vout(2, :)];  % polarization multiplexed signal    
            elseif size(Vin, 1) == 1 && not(isreal(Vin)) 
            %% Single-pol I/Q
                % Breaks data into in-phase, quadrature
                Vxi = real(Vin(1, :))/self.Vpi + self.Vbias;
                Vxq = imag(Vin(1, :))/self.Vpi + self.Vbias;

                % Filter driving signals
                Vxi = real(ifft(fft(Vxi).*Hmod));
                Vxq = real(ifft(fft(Vxq).*Hmod)); 

                % Modulate signal fields (each has unit average power)
                Enorm = 1/sqrt(sqrt(2)); % normalize so that E(|Vout|^2) = 1
                Vout   = Enorm*(sin(pi*Vxi/2) + 1i*sin(pi*Vxq/2));

                Eout =  Ein.*Vout(1, :);  % polarization multiplexed signal
            elseif size(Vin, 1) == 1 && isreal(Vin) 
            %% Single-pol intensity modulator
                % Breaks data into in-phase, quadrature in two pols
                Vx = Vin/self.Vpi + self.Vbias;

                % Filter driving signal
                Vx = real(ifft(fft(Vx).*Hmod));

                % Modulate signal fields (each has unit average power)
                Vout  = sin(pi*Vx/2);    

                Eout =  Ein.*Vout;  % polarization multiplexed signal
            else
                error('mzm: Invalid driving signal. Vin must be either 1 x N (real or complex) or 2 x N (complex)')
            end
        end
        
        function H = mismatch(self, f, verbose)
            %% Mach-Zehnder frequency response limited by velocity mismatch and loss
            % K. P. Ho, Phase-Modulated Optical Communication Systems. New York: Springer, 2005.
            % Inputs:
            % - ratio: velocity mismatch ratio e.g., 0.98
            % - L: Iteractive length in m
            % - f: frequency vector
            % - verbose (optional, default = false): whether to plot figure
            % Output:
            % - Mod: struct containing MZM parameters
                   
            omega = 2*pi*f;                    
            d_12   = self.n_r*(self.ratio-1)/self.c;                               % Velocity mismatch difference between optical and electrical waveguides
            a      = 0.01*100*sqrt(abs(omega)/2/pi*1e-9);                          % Microwave attenuation coefficient (electrode loss)
            H    = (1-exp(self.L*(-a + 1j*d_12*omega)))./(a-1j*d_12*omega)/self.L; % Freq. response of the optical modulator limited by elec. loss and velocity mismatch               
            H(isnan(H)) = 1;                                % Mod.Hel(f=0) is NaN 
             
            [H, delay] = remove_group_delay(H, f); % remove group delay

            if exist('verbose', 'var') && verbose
                figure(111), hold on
                plot(f/1e9, abs(Mod.H).^2)
                plot(f([1 end])/1e9, [0.5 0.5], 'k')
                plot(Mod.BW/1e9, 0.5, 'xk', 'MarkerSize', 6)
                xlabel('Frequency (GHz)')
                ylabel('|H(f)|^2')
                legend('MZM frequency response', sprintf('BW = %.2f GHz', Mod.BW/1e9));
                title('MZM frequency response')
            end
        end
                
    end
end
    