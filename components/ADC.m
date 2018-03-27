classdef ADC
    %% Analog-to-digital converter
    properties
        T % sampling period (in samples)
        ENOB = Inf; % effective number of bits
        Filter % instance of class Filter used as anti-aliasing filter
        offset = 0 % additional time offset (in number of samples) to be removed i.e., x(t) = x(t - offset/fs). Need not be integer
        excursion = [] % excursion limits. Specify as a [xmin, xmax]
        rclip = 0 % clipping ratio. Clipping ratio is defined as the percentage of the signal amplitude that is clipped in both extremes (see code)
    end
    
    methods
        function obj = ADC(T, ENOB, Filter)
            %% Constructor
            obj.T = T;
            obj.ENOB = ENOB;
            obj.Filter = Filter;
        end
        
        function [xq, xs, xa] = convert(self, x, f, timeRefSignal)
            %% Main method: Anti-aliasing filter -> sampling -> quantize
            xa = self.aa_filter(x, f);
            if exist('timeRefSignal', 'var') % check whether time reference was passed
                xs = self.sample(xa, timeRefSignal);
            else
                xs = self.sample(xa);
            end
            xq = self.quantize(xs);
        end
        
        function xa = aa_filter(self, x, f)
            %% Filter by antialiasing filter defined in class properties (Filter)
            % f is normalized frequency from -0.5 to 0.5            
            if self.offset ~= 0 % time shift in number of samples. Need not be integer
                Hshift = ifftshift(exp(1j*2*pi*f*self.offset));
            else
                Hshift = 1;
            end
            
            Haa = ifftshift(self.Filter.H(f));
            xa = real(ifft(fft(x).*Haa.*Hshift));           
        end
        
        function xs = sample(self, xa, timeRefSignal)
            %% Sample signal xa with period T (in samples). 
            % Inputs:
            % - xa: signal to be sampled (downsampled)
            % - T: sampling period in samples. T does not need to be an integer
            % - timeRefSignal (optional): reference signal. xa will be
            % aligned to timeRefSignal before downsampling
            
            % Align to time reference
            if exist('timeRefSignal', 'var')
                [c, lags] = xcorr(timeRefSignal, xa);
                [~, idx] = max(abs(c));
                if lags(idx) ~= 0
                    xa = circshift(xa, [0 lags(idx)]);
                    fprintf('ADC: input signal was delayed by %d samples to match time reference signal\n', lags(idx));
                end
            end

            % Downsample
            if isInteger(self.T) % resampling is not required
                xs = xa(1:self.T:end);
            else
                [N, D] = rat(self.T);
                fprintf('ADC: T is not interger, so signal was resampled by %d/%d.\n', N, D);
                xs = resample(xa, N, D);
            end
            
        end
        
        function [xq, varQ] = quantize(self, xs)
            %% Quantization and clipping
            if not(isinf(self.ENOB))
                if not(isempty(self.excursion))
                    xmax = self.excursion(2);
                    xmin = self.excursion(1);
                else
                    xmax = max(xs);
                    xmin = min(xs);
                end
                xamp = xmax - xmin;    

                % Clipping
                % clipped: [xmin, xmin + xamp*rclip) and (xmax - xamp*rclip, xmax]
                % not clipped: [xmin + xamp*rclip, xmax - xamp*rclip]
                xmin = xmin + xamp*self.rclip; % discounts portion to be clipped
                xmax = xmax - xamp*self.rclip;
                xamp = xmax - xmin;

                dx = xamp/(2^(self.ENOB)-1);
                Pc = mean(xs > xmax | xs < xmin);
                if Pc ~= 0
                    fprintf('ADC: clipping probability = %G\n', Pc) 
                end

                codebook = xmin:dx:xmax;
                partition = codebook(1:end-1) + dx/2;
                [~, xq, varQ] = quantiz(xs, partition, codebook); 
                fprintf('ADC: quantization noise variance = %G | Signal-to-quantization noise ratio = %.2f dB\n', varQ, 10*log10(var(xs)/varQ));
            else
                xq = xs;
                varQ = 0;
            end
        end
    end
end