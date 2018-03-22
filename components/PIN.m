%% Positive-intrinsic-negative photodiode
classdef PIN < APD % inherits methods and properties from apd
    %% Positive-negative photodiode inherits methods and properties from class apd
    properties
        % defined in apd.m
        % R    % responsivity
        % Id   % dark current  
        % BW   % bandwidth
    end
                
        %% Methods inherited from APD.m     
        % function Hapd = H(this, f) % APD frequency response. Normalized to have unit gain at DC
        % function hapd = ht(this, t) % APD impulse response       
        % function [Hw, y] = Hwhitening(self, f, P, N0, x) % Design whitening filter        
        % function bw = noisebw(this) % Noise bandwidth (one-sided)
        % function sig2 = varShot(this, Pin, Df) % Shot noise variance
        % function output = detect(this, Ein, fs, noise_stats, N0) % Direct detection              
        % function noise_std = stdNoise(this, Hrx, Hff, N0, RIN, sim) % Returns function handle noise_std, which calculates the noise standard deviation for a given power level P. 
        % function [Gopt, mpam] = optGain(this, mpam, tx, fiber, rx, sim) % Optimize APD gain: Given target BER finds APD gain that leads to minimum required optical power     
        % function [Pmean, mpam] = optimize_PAM_levels(this, Gapd, mpam, Tx, Fiber, Rx, sim) % Calculate optimal level spacing for a given APD gain     
        
        % Get and Set Methods defined in APD.m
        % function Fa = get.Fa(this) % excess noise factor      
        % function GaindB = get.GaindB(this) % Calculate gain in dB        
        % function BW = get.BW(this) % APD bandwidth
        % function Geff = get.Geff(this) % Effective gain = Responsivity x Gain
        % function this = set.GaindB(this, GdB) % Set gain in dB
        % function self = setGain(self, Gain) % Set APD Gain
        
    methods
        function self = PIN(R, Id, BW)
            %% Constructor
            % Inputs:
            % - R (optional, default = 1 A/W) = responsivity
            % - Id (optional defualt = 0 nA) = dark current (A)       
            % - BW (optional, default = Inf) = bandwidth of 1st-order
            % frequency response
     
            if ~exist('R', 'var')
                R = 1;
            end
            
            if ~exist('BW', 'var')
                BW = Inf;
            end
            
            if ~exist('Id', 'var')
                Id = 0;
            end
            % Call constructor from APD
            % apd(GaindB, ka, BW, R, Id)
            self@APD(0, 0, BW, R, Id);
        end       
    
        function PINtable = summary(self)
            %% Generate table summarizing class values
            disp('-- PIN class parameters summary:')
            rows = {'Responsivity'; 'Dark current'; 'Bandwidth'};
            Variables = {'R'; 'Id'; 'BW'};
            Values = [self.R; self.Id*1e9; self.BW/1e9];
            Units = {'A/W'; 'nA'; 'GHz'};

            PINtable = table(Variables, Values, Units, 'RowNames', rows);
        end
    end
end