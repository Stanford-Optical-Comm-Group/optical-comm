classdef Filter < handle
    %% Filter
    % Filter can be given in terms of num and den coefficients or frequency
    % response H(f) 
    % This class can be used as a filter within a loop. It can take one 
    % value at each iteration. This is useful for simulating a PLL. When
    % used this way, the filter is implemented as a direct form I. 
    % If a vector is passed as a parameter, then the filtering operation
    % will be performed using the filter function from Matlab (transpose
    % direct form I), or if the flag removeGroupDelay is set, the filter
    % will be implemented in the frequency domain
    properties 
        num % numerator = num(1) + num(2)z^-1 + num(3)z^-2+...
        den % denominator = den(1) + den(2)z^-2 + den(3)z^-2+...
        fs = 1 % sampling frequency 
        removeGroupDelay = true; % whether group delay should be removed. This is only used when a vector is passed to function filter
        % Note: group delay is calculated as the group delay at DC
    end
    
    properties(Dependent)
        filterType % FIR or IIR
        BW
    end
    
    properties
        memForward % register for forward path
        memFeedback % register for feedback path
    end
    
    properties(Access=private, Hidden)
        Hf = [];
    end
    
    methods
        function obj = Filter(varargin)
            %% Constructor  
            if nargin == 0 
                %% > Filter() % empty class
                obj.num = [];
                obj.den = [];  
            elseif not(isreal(varargin{1})) 
                %% > Filter(H, f) % implements filter given in H (num & den are unknown)
                obj.Hf = varargin{1};
                obj.num = [];
                obj.den = [];
            elseif ischar(varargin{1}) 
                %% > Filter(type, order, cutoff) % Filter specifications are passed: type, order, normalized cutoff frequency, fs (optional)
                % Example: Filter('butter', 5, 0.5)
                type = varargin{1};
                order = varargin{2};
                fcnorm = varargin{3};
                if length(varargin) == 4
                    obj.fs = varargin{4};
                else
                    obj.fs = 1;
                end
                
                filt = design_filter(type, order, fcnorm);
                obj.num = filt.num;
                obj.den = filt.den;
            elseif isstruct(varargin{1}) 
                %% > Filter(filt) % Filter struct designed with design_filter.m (legacy)
                filt = varargin{1};
                if nargin == 2
                    obj.fs = varargin{2};
                else
                    obj.fs = 1;
                end
                
                obj.num = filt.num;
                obj.den = filt.den;              
            else 
                %% > Filter(num, den) % Filter is already known: passing num, den, fs (optional)
                obj.num = varargin{1};
                obj.den = varargin{2};

                if nargin == 3
                    obj.fs = varargin{3};
                else
                    obj.fs = 1;
                end
            end
            
            obj.memForward = zeros(1, length(obj.num));
            obj.memFeedback = zeros(1, length(obj.den)-1);
        end
        
        function y = filter(self, x)
            %% Filter function
            if length(x) > 1 % x is the entire signal so can filter all at once
                if self.removeGroupDelay % whether to remove group delay
                    delay = grpdelay(self.num, self.den, 1);
                    if isInteger(delay) && isempty(self.Hf)
                        % If group delay is integer and num & den are known
                        % filter and remove group delay by dealying output
                        % signal
                        y = filter(self.num, self.den, x);
                        if size(y, 1) > size(y, 2)
                            y = circshift(y, [delay 0]);
                        else
                           y = circshift(y, [0 delay]);
                        end
                    else
                        % Otherwise, implements filter in frequency domain
                        f = freq_time(length(x), self.fs);
                        y = ifft(fft(x).*ifftshift(self.H(f)));
                        if isreal(x) % if input is real, remove spurious imaginary component
                            y = real(y);
                        end
                    end
                else % filters using transpose direct form 1
                    y = filter(self.num, self.den, x);
                end
                return
            end

            self.memForward = self.shiftState(self.memForward, x);
            if strcmpi(self.filterType, 'FIR')
                y = sum(self.memForward.*self.num);
            else
                y = (sum(self.memForward.*self.num)- sum(self.memFeedback.*self.den(2:end)))/self.den(1);
                self.memFeedback = self.shiftState(self.memFeedback, y);
            end
        end
        
        function groupDelay = groupDelay(self)
            %% Override group delay method from parent class because only one filtering operation is performed
            if not(isempty(self.Hf)) % filter is given by Hf frequency response
                N = length(self.Hf);
                df = self.fs/N;
                f = -self.fs/2:df:self.fs/2-df;
                [~, groupDelay] = remove_group_delay(self.Hf, f);
            else % filter specified by coefficients num and den 
                groupDelay = grpdelay(self.num, self.den, 1)/self.fs;
            end
        end
        
        function fvtool(self)
            %% Opens Matlab's fvtool for this filter
            fvtool(self.num, self.den)
        end
        
        function Hf = H(self, f)
            %% Frequency response
            if isempty(self.Hf) % filter specified by coefficients num and den 
                Hf = freqz(self.num, self.den, f, self.fs);
                if self.removeGroupDelay
                    Hf = Hf.*exp(1j*2*pi*f*self.groupDelay);
                end
            else % filter specified by frequency response Hf
                if self.removeGroupDelay
                    N = length(self.Hf);
                    df = self.fs/N;
                    f = -self.fs/2:df:self.fs/2-df;
                    Hf = remove_group_delay(self.Hf, f);
                else
                    Hf = self.Hf;
                end
            end
        end
        
        function self = reset(self)
            %% Reset all states
            self.memForward = zeros(size(self.memForward));
            self.memFeedback = zeros(size(self.memFeedback));
        end
        
        function newFilter = cascade(this, otherFilter)
            %% Cascade this filter with another filter
            newNum = conv(this.num, otherFilter.num);
            newDen = conv(this.den, otherFilter.den);
            newFilter = Filter(newNum, newDen, this.fs);
        end
        
        function state = shiftState(~, state, x)
            %% Shift states for a new input
            state(2:end) = state(1:end-1);
            state(1) = x;
        end       
        
        function outputSignal = validate(self, inputSignal)
            %% Check if filter implementation is working properly
            if not(exist('inputSignal', 'var'))
                inputSignal = randn(1, 1000); % generates random signal for testing
            end
            
            self.reset();
            outputSignal = zeros(size(inputSignal));
            self = self.reset();
            for t = 1:length(inputSignal)
                outputSignal(t) = self.filter(inputSignal(t));
            end
            
            figure, hold on, box on
            plot(inputSignal)
            plot(outputSignal, 'r')
            plot(filter(self.num, self.den, inputSignal), '--k');
            legend('Input signal', 'Output signal', 'Output signal using matlab function filter')
            self.reset();
        end  
    end
    
    %% Get and set methods
    methods
        function filterType = get.filterType(self)
            if isempty(self.den) && isempty(self.num)
                filterType = 'unknown';
            elseif length(self.den) == 1
                filterType = 'FIR';
            else
                filterType = 'IIR';
            end
        end
    end
end