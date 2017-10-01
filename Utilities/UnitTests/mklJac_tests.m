%% A Unit Test Class for mklJac
classdef mklJac_tests < matlab.unittest.TestCase

    properties
        absTol = 1e-8;
    end
    
    % Unit Tests
    methods (Test)
        
        %-- Input Args --%
        function inputArgs(testCase)
            testCase.verifyError(@() mklJac(@(x) sin(x)), 'OPTIMex:InputError');   % not enough args
            testCase.verifyError(@() mklJac(1, 1), 'OPTIMex:InputError');   % not fcn handle
            testCase.verifyError(@() mklJac(@(x) sin(x), int16(1)), 'OPTIMex:InputError'); % wrong input type
            testCase.verifyError(@() mklJac(@(x) sin(x), 1i), 'OPTIMex:InputError'); % wrong input type
            testCase.verifyError(@() mklJac(@(x) sin(x), [1 1; 1 1]), 'OPTIMex:InputError'); % wrong input type
            testCase.verifyError(@() mklJac(@(x) sin(x), 1, int16(1)), 'OPTIMex:InputError'); % wrong input type
            testCase.verifyError(@() mklJac(@(x) sin(x), 1, [1;1]), 'OPTIMex:InputError'); % wrong input type
            testCase.verifyError(@() mklJac(@(x) sin(x), 1, 0), 'OPTIMex:InputError'); % wrong val
            testCase.verifyError(@() mklJac(@(x) sin(x), 1, 1e9), 'OPTIMex:InputError'); % wrong val
            testCase.verifyError(@() mklJac(@(x) sin(x), 1, 1, int16(1)), 'OPTIMex:InputError'); % wrong input type
            testCase.verifyError(@() mklJac(@(x) sin(x), 1, 1, [1;1]), 'OPTIMex:InputError'); % wrong input type
            testCase.verifyError(@() mklJac(@(x) sin(x), 1, 1, 1e-18), 'OPTIMex:InputError'); % wrong val
            testCase.verifyError(@() mklJac(@(x) sin(x), 1, 1, 1.1), 'OPTIMex:InputError'); % wrong val
            testCase.verifyError(@() mklJac(@(x) sin(x), [1;1], 1), 'OPTIMex:InputError'); % wrong length
            testCase.verifyError(@() mklJac(@(x) sin(x), [1;1], 3), 'OPTIMex:InputError'); % wrong length
            testCase.verifyError(@() mklJac(@(x) sin(x), 1, 2), 'OPTIMex:InputError'); % wrong length
        end
        
        %-- Valid Operation --%
        function scalarDiff(testCase)
            testCase.verifyEqual(1, mklJac(@(x) sin(x), 0), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(-1, mklJac(@(x) sin(x), pi), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(0, mklJac(@(x) cos(x), pi), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(0, mklJac(@(x) cos(x), -pi), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(-1, mklJac(@(x) cos(x), pi/2), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(1, mklJac(@(x) cos(x), -pi/2), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(2, mklJac(@(x) 2*x, 0), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(2, mklJac(@(x) 2*x, 1), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(2, mklJac(@(x) 2*x, 2), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(-2, mklJac(@(x) -2*x, 0), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(-2, mklJac(@(x) -2*x, 1), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(-2, mklJac(@(x) -2*x, 2), 'AbsTol', testCase.absTol);
        end
        
        function scalarDiffSweep(testCase)
            fun = @(x) 3*cos(x^2) + 0.5*sin(x/2);
            grad = @(x) cos(x/2)/4 - 6*x*sin(x^2);
            x = linspace(-pi,pi);
            for i = 1:length(x)
                testCase.verifyEqual(grad(x(i)), mklJac(fun, x(i)), 'AbsTol', testCase.absTol);
            end
        end
        
        function vectorDiff(testCase)
            testCase.verifyEqual(eye(3), mklJac(@(x) sin(x), zeros(3,1)), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(eye(6), mklJac(@(x) sin(x), zeros(6,1)), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(-eye(3), mklJac(@(x) sin(x), pi*ones(3,1)), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(zeros(3), mklJac(@(x) sin(x), pi/2*ones(3,1)), 'AbsTol', testCase.absTol);
        end
        
        function vectorDiffSweep(testCase)
            fun = @(x) 3*cos(x(1)^2)*sin(x(2)) + 0.5*sin(x(2)/2);
            grad = @(x) [ -6*x(1)*sin(x(1)^2)*sin(x(2)), cos(x(2)/2)/4 + 3*cos(x(1)^2)*cos(x(2))];
            x1 = linspace(-pi,pi);
            x2 = linspace(pi,-pi);
            for i = 1:length(x1)
                testCase.verifyEqual(grad([x1(i);x2(i)]), mklJac(fun, [x1(i);x2(i)]), 'AbsTol', testCase.absTol);
            end
        end
        
        function vectorFunDiffSweep(testCase)
            fun = @(x) [100*(x(2)-x(1)^2); 1 - x(1)];
            grad = @(x) [-200*x(1) 100; -1 0];
            x = linspace(-1,1,10);
            for i = 1:length(x)
                for j = 1:length(x)
                    testCase.verifyEqual(grad([x(i);x(j)]), mklJac(fun, [x(i);x(j)]), 'AbsTol', testCase.absTol*10); %this is a hard one
                end
            end
        end
        
        function autoSizeIdentify(testCase)
            testCase.verifyEqual(1, numel(mklJac(@(x) sin(x), zeros(1,1))));
            testCase.verifyEqual(4, numel(mklJac(@(x) sin(x), zeros(2,1))));
            testCase.verifyEqual(9, numel(mklJac(@(x) sin(x), zeros(3,1))));
            testCase.verifyEqual(16, numel(mklJac(@(x) sin(x), zeros(4,1))));
        end                
    end
    
end

