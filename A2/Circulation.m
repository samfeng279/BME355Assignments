classdef Circulation < handle
    % Model of systemic circulation from Ferreira et al. (2005), A Nonlinear 
    % State-Space Model of a Combined Cardiovascular System and a Rotary Pump, 
    % IEEE Conference on Decision and Control. 
    % 
    % NOTE: Several parts are left incomplete. Search the text for "TODO" to find them.  
    
    properties (SetAccess = private)
        % These shouldn't be set individually because their values are
        % related, so we have a setter function that changes them all
        % together. 
        HR;
        tc;
        Tmax;
    end
    
    properties (Access = public)
        Emax;
        Emin;
        
        nonSlackBloodVolume = 210; % ml
        
        R1 = .5; % between .5 and 2
        R2 = .005;
        R3 = .001;
        R4 = .0398;

        C2 = 4.4;
        C3 = 1.33;

        L = .0005;
    end
    
    methods (Access = public)
        function C = Circulation(HR, Emax, Emin) 
            % HR: heart rate (beats per minute)
            % Emax: maximum elastance
            % Emin: minimum elastance
            
            C.setHeartRate(HR);
            C.Emax = Emax;
            C.Emin = Emin;
        end
        
        function setHeartRate(C, HR)
            % HR: heart rate (beats per minute)
            
            % Recall C (first argument) is the Circulation instance. 
            
            C.HR = HR;
            C.tc = 60/HR; 
            C.Tmax = .2+.15*C.tc; % contraction time
        end
        
        function dx = getDerivative(C, t, x)
            % t: time
            % x: state variables [ventricular pressure; atrial pressure; arterial pressure; aortic flow]
            % dx: time derivatives of state variables
            
            flow = x(4);

            if x(2) > x(1)
                A = filling(C, t);
            elseif % TODO: fill in this expression 
                A = ejection(C, t);
            else 
                A = isovolumic(C, t);
            end

            dx = A*x;
        end
        
        function A = isovolumic(C, t)
            % This method produces the isovolumic A matrix from the paper.
            % 
            % t: time (needed because elastance is a function of time)
            
            el = elastance(C,t);
            del_dt = elastanceFiniteDifference(C, t);
            A = [del_dt/el 0 0 0;
                0 -1/(C.R1*C.C2) 1/(C.R1*C.C2) 0; 
                0 1/(C.R1*C.C3) -1/(C.R1*C.C3) 0; 
                0 0 0 0];
        end
        
        function A = ejection(C, t)
            % TODO: implement this method 
        end
        
        function A = filling(C, t)
            % TODO: implement this method 
        end
        
        function result = elastance(C, t)
            % t: time (needed because elastance is a function of time)  
            % result: time-varying elastance
            
            % C.Tmax is the length of time of the ventricular contraction 
            % tn is time normalized to C.Tmax
            tn = rem(t,C.tc)/C.Tmax;
            neg = find(tn < 0); %this line and the next are for generality 
            tn(neg) = tn(neg) + C.Tmax; % in case you give a t < 0
            En = 1.55 * (tn/.7).^1.9 ./ (1 + (tn/.7).^1.9) ./ (1 + (tn/1.17).^21.9); %En is normalized elastance (between 0 and 1)
            result = (C.Emax-C.Emin)*En+C.Emin; % this rescales the normalized elastance to the range (Emin,Emax)
        end

        function result = elastanceDerivative(C, t)
            % t: time (needed because elastance is a function of time)  
            % result: time derivative of time-varying elastance            
            
            tn = rem(t,C.tc)/C.Tmax;
            if tn < 0, tn = tn + C.Tmax; end;
            
            num = 1.55 * (tn/.7).^1.9;
            den = (1 + (tn/.7).^1.9) .* (1 + (tn/1.17).^21.9);
                        
            dnum_dtn = 1.55 * 1.9 * (tn/.7).^0.9 * (1/.7);
            dden_dtn = (21.9/1.17)*(tn/1.17).^20.9 .* (1 + (tn/.7).^1.9) + (1.9/.7)*(tn/.7).^0.9 .* (1 + (tn/1.17).^21.9);
            
            dtn_dt = 1/C.Tmax;
            % this is the quotient rule applied to numerical values of the 
            % numerator, denominator, and their derivatives
            dEn_dtn = ((den .* dnum_dtn) - (num .* dden_dtn)) ./ den.^2; 
            dE_dEn = (C.Emax-C.Emin);
            result = dE_dEn * dEn_dtn * dtn_dt;
        end
        
        function result = elastanceFiniteDifference(C, t)
            % t: time (needed because elastance is a function of time)  
            % result: finite-difference approximation of time derivative of 
            %   time-varying elastance            

            dt = .0001;
            forwardTime = t + dt;
            backwardTime = max(0, t - dt); % small negative times are wrapped to end of cycle
            forward = elastance(C, forwardTime);
            backward = elastance(C, backwardTime);
            result = (forward - backward) ./ (forwardTime - backwardTime);
        end
        
        function checkDerivative(C)
            % This method checks the finite difference against the symbolic
            % derivative and plots both if they are very different. 
            
            t = 0:.005:1; 
            d = elastanceDerivative(C, t);
            fd = elastanceFiniteDifference(C, t);
            if max(abs(d - fd)) > .1          
                plot(t, d, 'r', t, fd, 'g--')
                title(sprintf('Largest difference: %f', max(abs(d - fd))));
            end
        end
        
        function [time, state] = simulate(C, totalTime)
            % totalTime: # seconds to simulate
            % time: times at which the state is estimated
            % state: state vector at each time
            
            % We put all the blood pressure in the atria as an
            % initial condition. Note that we can't get the total blood
            % volume by multiplying by C2, because we're missing the
            % pulmonary loop. 
            
            initialState = [0; C.nonSlackBloodVolume/C.C2; 0; 0]; 
            ofun = @(t,x) C.getDerivative(t,x); % wrapper function because ode45 expects a function rather than a method
            
            % TODO: get time and state with a call to ode45
        end
        
    end
        
end