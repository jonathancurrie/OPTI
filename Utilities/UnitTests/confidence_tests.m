%% A Unit Test Class for confidence interval and fit statistics
classdef confidence_tests < matlab.unittest.TestCase

    properties
        absTol          = 1e-10;
        covAbsTol       = 1e-6;
        confIntAbsTol   = 1e-6;
        confBndsRelTol  = 1e-6;
        pValAbsTol      = 1e-8;
        stdErrRelTol    = 1e-5;
        tStatRelTol     = 1e-5;
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
        
        %-- Himmelblau - with weights --%
        function himmelblauWeights(testCase)
            mdl = testCase.data.himmelw;
            mdl.fcn = @confidence_tests.himmelblauFcn;
            % Check fit statistics
            checkFitStats(testCase,mdl);
        end
        
        %-- Himmelblau - with lower bound active --%
        function himmelblauLb(testCase)
            mdl = testCase.data.himmellb;
            mdl.fcn = @confidence_tests.himmelblauFcn;
            % Check fit statistics
            checkFitStats(testCase,mdl);
        end
        
        %-- Himmelblau - with upper bound active --%
        function himmelblauUb(testCase)
            mdl = testCase.data.himmelub;
            mdl.fcn = @confidence_tests.himmelblauFcn;
            % Check fit statistics
            checkFitStats(testCase,mdl);
        end
        
        %-- Himmelblau - with both bounds active --%
        function himmelblauBothBounds(testCase)
            mdl = testCase.data.himmelbnd;
            mdl.fcn = @confidence_tests.himmelblauFcn;
            % Check fit statistics
            checkFitStats(testCase,mdl);
        end
        
        %-- Himmelblau - with both bounds and weights --%
        function himmelblauBothBoundsWeights(testCase)
            mdl = testCase.data.himmelbndw;
            mdl.fcn = @confidence_tests.himmelblauFcn;
            % Check fit statistics
            checkFitStats(testCase,mdl);
        end
        
        %-- SAS - no weights --%
        function sas(testCase)
            mdl = testCase.data.sas;
            mdl.fcn = @confidence_tests.sasFcn;
            % Check fit statistics
            checkFitStats(testCase,mdl);
        end
        
        %-- SAS - with weights --%
        function sasWeights(testCase)
            mdl = testCase.data.sasw;
            mdl.fcn = @confidence_tests.sasFcn;
            % Check fit statistics
            checkFitStats(testCase,mdl);
        end
    end
    
    methods
       function checkFitStats(testCase, mdl)
            % Build OPTI object
            optiObj = opti('fun',mdl.fcn, 'data', mdl.xdata, mdl.ydata, 'weights', mdl.wts, 'x0', mdl.sol, 'bounds', mdl.lb, mdl.ub);
            % Manually set the solution
            optiObj.setSolution(mdl.SSE, mdl.sol);
            % Calculate the fit statistics;
            stats  = calcStatistics(optiObj, 0.95, false);
            
            % Check statistcs
            testCase.verifyEqual(mdl.Rsquare, stats.Rsquare, 'AbsTol', testCase.absTol);
            testCase.verifyEqual(mdl.AdjRsquare, stats.AdjRsquare, 'AbsTol', testCase.absTol);
            testCase.verifyEqual(mdl.RMSE, stats.RMSE, 'AbsTol', testCase.absTol);
            testCase.verifyEqual(mdl.DFE, stats.DFE, 'AbsTol');
            if (isfield(mdl, 'cov'))
                testCase.verifyEqual(mdl.cov, stats.Cov, 'AbsTol', testCase.covAbsTol);
            end
            testCase.verifyEqual(mdl.confInt, stats.ConfInt, 'AbsTol', testCase.confIntAbsTol);
            testCase.verifyEqual(mdl.confBnds, stats.ConfBnds.bnds, 'RelTol', testCase.confBndsRelTol);
            if (isfield(mdl, 'param'))
                testCase.verifyEqual(mdl.param(:,2), stats.Param.StdError, 'RelTol', testCase.stdErrRelTol);
                testCase.verifyEqual(mdl.param(:,3), stats.Param.tStat, 'RelTol', testCase.tStatRelTol);
                testCase.verifyEqual(mdl.param(:,4), stats.Param.pValues, 'AbsTol', testCase.pValAbsTol);
            end
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

