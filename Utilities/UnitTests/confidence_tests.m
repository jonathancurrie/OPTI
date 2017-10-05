%% A Unit Test Class for confidence interval and fit statistics
classdef confidence_tests < matlab.unittest.TestCase

    properties
        absTol = 1e-10;
        data   = load('Utilities/UnitTests/UnitTestData/confidenceData.mat');
    end
    
    % Unit Tests
    methods (Test)

        %-- Himmelblau - no weights --%
        function himmelblau(testCase)
            mdl = testCase.data.himmel;
            mdl.fcn = @confidence_tests.himmelblauFcn;
            % Check fit statistics
            checkFitStats(testCase,mdl);
        end
    end
    
    methods
       function checkFitStats(testCase, mdl)
            % Build OPTI object
            optiObj = opti('fun',mdl.fcn, 'data', mdl.xdata, mdl.ydata, 'weights', mdl.wts, 'x0', mdl.param(:,1));
            sol = mdl.param(:,1);
            fval = sum((mdl.fcn(sol,mdl.xdata) - mdl.ydata).^2);
            optiObj.setSolution(fval, sol);
            stats  = optiObj.calcStatistics
            
            testCase.verifyEqual(f.qf, rmathlib('qf',f.p,f.df1,f.df2), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(f.df, rmathlib('df',f.x,f.df1,f.df2), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(f.pf, rmathlib('pf',f.q,f.df1,f.df2), 'AbsTol', testCase.absTol);
        end
    end
    
    
    methods (Static)
        function r = himmelblauFcn(theta,p)
            % r = f(x,theta) (x = pressure, theta unknowns)
            r = theta(1)*p./(1+theta(2)*p); 
        end
        
        function k = sasFcn(theta,n)
            k = theta(1)*n./(theta(2) + n);
        end        
    end
end

