classdef DAC
    %% Digital-to-analog conversion
    properties
        Nhold % upsampling factor (integer)
        resolution % resolution in bits
        Filter % instance of class filter corresponding to DAC's analog frequency response
        offset = 0 % additional time offset (in number of samples) to be removed. Need not be integer
        excursion = [] % excursion limits. Specify as a [xmin, xmax]
        rclip = 0 % clipping ratio. Clipping ratio is defined as the percentage of the signal amplitude that is clipped in both extremes (see code)
    end
    
    methods
        function obj = DAC(Nhold, resolution, Filter)
            obj.Nhold = Nhold;
            obj.resolution = resolution;
            obj.Filter = Filter;
        end
        
        function xa = convert(self, x, f)
            %% Main method
            xq = self.quantize(x);
            xa = self.sample_and_hold(xq, f);
        end
        
        function [xq, varQ] = quantize(self, x)
            %% Quantization and clipping
            if not(isinf(self.resolution))
                if not(isempty(self.excursion))
                    xmax = self.excursion(2);
                    xmin = self.excursion(1);
                else
                    xmax = max(x);
                    xmin = min(x);
                end
                xamp = xmax - xmin;    

                % Clipping
                % clipped: [xmin, xmin + xamp*rclip) and (xmax - xamp*rclip, xmax]
                % not clipped: [xmin + xamp*rclip, xmax - xamp*rclip]
                xmin = xmin + xamp*self.rclip; % discounts portion to be clipped
                xmax = xmax - xamp*self.rclip;
                xamp = xmax - xmin;

                dx = xamp/(2^(self.resolution)-1);
                Pc = mean(x > xmax | x < xmin);
                if Pc ~= 0
                    fprintf('DAC: clipping probability = %G\n', Pc) 
                end

                codebook = xmin:dx:xmax;
                partition = codebook(1:end-1) + dx/2;
                [~, xq, varQ] = quantiz(x, partition, codebook); 
                fprintf('DAC: quantization noise variance = %G | Signal-to-quantization noise ratio = %.2f dB\n', varQ, 10*log10(var(x)/varQ));
            else
                xq = x;
                varQ = 0;
            end
        end
        
        function xa = sample_and_hold(self, xq, f)
            %% Zero-order holder (ZOH)
            assert(isInteger(self.Nhold), 'DAC: oversampling ratio of DAC must be an integer.');
            xqzoh = upsample(xq, self.Nhold);
            xqzoh = filter(ones(1, self.Nhold), 1, xqzoh); % group delay due to ZOH is removed inn the frequency domain

            % Filter
            Hdac = ifftshift(self.Filter.H(f)); % filter frequency response

            % Time shift. Removes any group delay specified by offset. This is
            % equivalent to resampling the signal
            nd = self.offset + (self.Nhold-1)/2; % if offtset==0, only removed offset due to ZOH

            Hshift = ifftshift(exp(1j*2*pi*f*nd));

            % Filtering
            xa = real(ifft(fft(xqzoh).*Hdac.*Hshift)); % filter
        end 
    end
end
    
            
