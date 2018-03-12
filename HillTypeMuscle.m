classdef HillTypeMuscle < handle
    % Damped Hill-type muscle model adapted from Millard et al. (2013). The
    % dynamic model is defined in terms of normalized length and velocity. 
    % To model a particular muscle, scale factors are needed for force, CE
    % length, and SE length. These are given as constructor arguments. 
    
    properties (Access = public)
        f0M;               
        restingLengthCE;    
        restingLengthSE;    
    end
    
    properties (Constant)
        % curve fits for CE force-length and force-velocity 
        forceLengthRegression = getCEForceLengthRegression();
        forceVelocityRegression = getCEForceVelocityRegression();
    end
    
    methods (Access = public)
        function m = HillTypeMuscle(f0M, restingLengthCE, restingLengthSE)
            % The arguments are scale factors model is normalized 
            % f0M: maximum isometric force
            % restingLengthCE: actual length of CE (m) that corresponds to 
            %   normalized length of 1
            % restingLengthSE: % actual length of SE (m) that corresponds 
            %   to normalized length of 1
            
            m.f0M = f0M;
            m.restingLengthCE = restingLengthCE;
            m.restingLengthSE = restingLengthSE;
        end
        
        function result = getNormalizedLengthSE(m, muscleTendonLength, normalizedLengthCE)
            % Calculates the normalized length of the series elastic
            % element. 
            % 
            % muscleTendonLength: physical (non-normalized) length of the
            %   full muscle-tendon complex (typically found from joint 
            %   angles and musculoskeletal geometry)
            % normalizedLengthCE: normalized length of the contractile
            %   element (the state variable of the muscle model)
            % result: normalized length of the series elastic element
            
            result = (muscleTendonLength - m.restingLengthCE*normalizedLengthCE) / m.restingLengthSE;
        end
        
        function result = getForce(m, length, normalizedLengthCE)
            % length: muscle-tendon length (m)
            % normalizedLengthCE: normalized length of CE (the state
            %   variable)
            result = m.f0M * m.forceLengthSE(m.getNormalizedLengthSE(length, normalizedLengthCE));
        end
                
        function simulateIsometric(m)
            L = m.restingLengthCE + m.restingLengthSE;
            afun = @(t) t>0.5;
            
            % WRITE CODE HERE (to define odefun)
            odefun = @(t, x) HillTypeMuscle.getVelocity(afun(t), x, m.getNormalizedLengthSE(L, x));
            
            OPTIONS = odeset('AbsTol', 1e-6, 'RelTol', 1e-5);  % use this as the final argument to ode45 
            [time, x] = ode45(odefun, [0 2], 1, OPTIONS);
            
            figure
            subplot(2,1,1)
            plot(time, x*m.restingLengthCE)
            ylabel('CE length')
            set(gca, 'FontSize', 18)
            subplot(2,1,2)
            plot(time, m.f0M*HillTypeMuscle.forceLengthSE(m.getNormalizedLengthSE(L, x)));
            xlabel('Time (s)')
            ylabel('Force')
            set(gca, 'FontSize', 18)
        end
        
    end
    
    methods (Static)
        
        function result = getVelocity(a, lM, lT)
            % Calculates normalized velocity of contractile element. 
            % 
            % a: activation (between 0 and 1)
            % lM: normalized length of contractile element (this is
            %   \tilde{l}^M in the paper and the lecture slides)
            % lT: normalized length of series elastic element (this is
            %   \tilde{l}^T in the paper and the lecture slides)
            % result: normalized velocity
            
            beta = 0.1; % damping coefficient (see damped model in Millard et al.)
            
            seFL = HillTypeMuscle.forceLengthSE(lT);
            ceFL = HillTypeMuscle.forceLengthCE(lM);
            peFL = HillTypeMuscle.forceLengthPE(lM);
            
            fun = @(vM) a*ceFL*HillTypeMuscle.forceVelocityCE(vM) + peFL + beta*vM - seFL;
            result = fzero(fun, 0); 
        end
        
        function result = forceLengthCE(lM)
            % Normalized force-length curve of contractile element. 
            % 
            % lM: contracile element length
            % result: force-length scale factor
            
            result = HillTypeMuscle.forceLengthRegression.eval(lM);            
        end
        
        function result = forceVelocityCE(vM)
            % Normalized force-velocity curve of contractile element.  
            % 
            % vM: contracile element velocity
            % result: force-velocity scale factor
            
            result = max(0, HillTypeMuscle.forceVelocityRegression.eval(vM));
        end
        
        function result = forceLengthSE(lT)
            % Normalized force-length curve of series elastic element. 
            % 
            % lT: normalized length of tendon (series elastic element)
            % result: force produced by tendon
            
            % WRITE CODE HERE
            sT = 1;
            result = zeros(size(lT));
            for i = 1:size(lT,2)
                if (lT(i) >= sT) 
                    result(i) = 10 * (lT(i) - sT) + 240 * (lT(i) - sT) ^ 2;
                else
                    result(i) = 0;
                end
            end
        end
        
        function result = forceLengthPE(lM)
            % Normalized force-length curve of parallel elastic element.
            % 
            % lM: normalized length of contractile element
            % result: force produced by parallel elastic element
            
            % WRITE CODE HERE            
            sM = 1;
            
            result = zeros(size(lM));
            for i = 1:size(lM,2)
                if (lM(i) >= sM) 
                    result(i) = 3 * (lM(i) - sM) ^ 2 / (0.6 + lM(i) - sM);
                else
                    result(i) = 0;
                end
            end
        end
        
        function plotCurves()
            % Plot force-length, force-velocity, SE, and PE curves. 
            
            lM = 0:.01:1.8;
            vM = -1.2:.01:1.2;
            lT = 0:.01:1.07;
            figure
            subplot(2,1,1), hold on
            plot(lM, HillTypeMuscle.forceLengthCE(lM), 'r')
            plot(lM, HillTypeMuscle.forceLengthPE(lM), 'g')
            plot(lT, HillTypeMuscle.forceLengthSE(lT), 'b')
            legend('CE', 'PE', 'SE', 'location', 'northwest')
            xlabel('Normalized length')
            ylabel('Force scale factor')
            set(gca, 'FontSize', 18)
            subplot(2,1,2)
            plot(vM, HillTypeMuscle.forceVelocityCE(vM), 'k')
            set(gca, 'FontSize', 18)
            xlabel('Normalized CE velocity')
            ylabel('Force scale factor')
        end
        
    end
    
end

function result = getCEForceVelocityRegression()
    % result: regression model of contractile element force-length curve
    %   based on data in Millard et al.
    
    data = [-1.0028395556708567, 0.0024834319945283845
    -0.8858611825192801, 0.03218792009622429
    -0.5176245843258415, 0.15771090304473967
    -0.5232565269687035, 0.16930496922242444
    -0.29749770052593094, 0.2899790099290114
    -0.2828848376217543, 0.3545364496120378
    -0.1801231103040022, 0.3892195938775034
    -0.08494610976156225, 0.5927831890757294
    -0.10185137142991896, 0.6259097662790973
    -0.0326643239546236, 0.7682365981934388
    -0.020787245583830716, 0.8526638522676352
    0.0028442725407418212, 0.9999952831301149
    0.014617579774061973, 1.0662107025777694
    0.04058866536166583, 1.124136223202283
    0.026390887007381902, 1.132426122025424
    0.021070257776939272, 1.1986556920827338
    0.05844673474682183, 1.2582274002971627
    0.09900238201929201, 1.3757434966156459
    0.1020023112662436, 1.4022310794556732
    0.10055894908138963, 1.1489210160137733
    0.1946227683309354, 1.1571212943090965
    0.3313459588217258, 1.152041225442796
    0.5510200231126625, 1.204839508502158];

    velocity = data(:,1);
    force = data(:,2);

    centres = linspace(-.5, 0, 6);
    width = .3;
    lambda = .01;
    result = Regression(velocity, force, centres, width, lambda, 1);
end

function result = getCEForceLengthRegression()
    % result: Regression model of normalized CE force-length curve. Data 
    %   from Winters et al. (2011) Fig. 3C, normalized so that max force is
    %   ~1 and length at max force is ~1. 
    
    % WRITE CODE HERE (use WebPlotDigitizer to extract force-length points 
    % from Winters et al. (2011) Figure 3C. Click "View Data", select all, 
    % cut, and paste below, then write code to normalize and perform the 
    % regression curve fit. Use the plotCurves() function to verify that 
    % the result is reasonable. 
    % 
    % Don't forget to normalize data so optimal length = 1 and peak = 1.
    
    data = [
        % start by pasting WebPlotDigitizer results here
        % sample active data
        40.34655403660477, 36.87001724439568
        45.71414290355635, 68.30080223927234
        56.87268304467716, 100.74225876590111
        57.79604961720776, 92.32499437682753
        53.347995234882006, 89.85579686601858
        67.38709919276235, 52.602904056181615
        72.84599171935788, 26.485392247519542
        75.28253317671758, 13.079749081548357
        62.63264439057306, 80.55581935870845
        43.76827531052409, 57.32736860520333
        40.31856313364823, 17.72423962212281
        37.265555361174286, 9.466923249943704
        50.777663925890764, 91.7491815160073
        54.65040528494907, 100.70427111188866
        65.71347645348597, 67.84495039112284
        42.09431934621248, 32.341489016069715
    ];
    
    [optForce, index] = max(data(:,2));
    optLength = data(index,1);
    normalizer = [optLength, optForce];
    data = bsxfun(@rdivide, data, normalizer);
    
    length = data(:,1);
    force = data(:,2);

    centres = linspace(0.6, 1.2, 6);
    width = .3;
    lambda = .01;
    result = Regression(length, force, centres, width, lambda, 1);
end
