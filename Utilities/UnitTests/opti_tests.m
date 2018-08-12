%% A Unit Test Class for the OPTI Object
classdef opti_tests < matlab.unittest.TestCase

    properties

    end
    
    % Unit Tests
    methods (Test)

        %-- Single Precision: Nonlinear --%
        function singleNonlinear(testCase)
            
            optiObj = opti_tests.makeHS71();
            defaultOpts = optiset(optiObj.opts, 'warnings', 'none', 'derivCheck', 'on');
            
            testCase.verifyError(@() opti(optiObj, 'fun', @(x) single(opti_tests.hs71_obj(x))), 'OPTI:NotDouble');
            testCase.verifyError(@() opti(optiObj, 'nlcon', @(x) single(opti_tests.hs71_con(x))), 'OPTI:NotDouble');           
            testCase.verifyError(@() opti(optiObj, 'grad', @(x) single(opti_tests.hs71_grad(x)), 'opts', defaultOpts), 'OPTI:NotDouble');         
            testCase.verifyError(@() opti(optiObj, 'nljac', @(x) single(opti_tests.hs71_jac(x)), 'opts', defaultOpts), 'OPTI:NotDouble');
            testCase.verifyError(@() opti(optiObj, 'hess', @(x,s,l) single(opti_tests.hs71_hess(x,s,l)), 'opts', defaultOpts), 'OPTI:NotDouble');
        end
    end
    
    methods (Static)
        %-- Hock & Schittkowski #71 --%
        function optiObj = makeHS71()            
            cl = [25;40];
            cu = [Inf;40];
            lb = ones(4,1);
            ub = 5*ones(4,1);        
            x0 = [1 5 5 1]';               
            opts = optiset('warnings', 'none');
            optiObj = opti('fun',@opti_tests.hs71_obj,'grad',@opti_tests.hs71_grad,'hess',@opti_tests.hs71_hess,...
                            'nl',@opti_tests.hs71_con,cl,cu,'nljac',@opti_tests.hs71_jac,'nljacstr',@opti_tests.hs71_jacStr,...
                            'hessstr',@opti_tests.hs71_hessStr,'bounds',lb,ub, 'x0', x0, 'opts', opts);
        end
            
        function obj = hs71_obj(x)
            obj = x(1)*x(4)*sum(x(1:3)) + x(3);
        end
        
        function grad = hs71_grad(x)
            grad = [ x(1)*x(4) + x(4)*sum(x(1:3))
                    x(1)*x(4)
                    x(1)*x(4) + 1
                    x(1)*sum(x(1:3)) ]';
        end
        
        function con = hs71_con(x)
            con = [ prod(x); sum(x.^2) ];
        end
        
        function jac = hs71_jac(x)
            jac = [ prod(x)./x'; 2*x' ];
        end
        
        function jacStr = hs71_jacStr()
            jacStr = sparse(ones(2,4));
        end
        
        function hess = hs71_hess(x, sigma, lambda)
            hess = sigma*[ 2*x(4)             0      0   0;
                              x(4)               0      0   0;
                              x(4)               0      0   0;
                              2*x(1)+x(2)+x(3)  x(1)  x(1)  0 ] + ...
                   lambda(1)*[   0          0         0         0;
                              x(3)*x(4)     0         0         0;
                              x(2)*x(4) x(1)*x(4)     0         0;
                              x(2)*x(3) x(1)*x(3) x(1)*x(2)     0  ] + ...
         		   lambda(2)*diag([2 2 2 2]);
        end
            
        function hessStr = hs71_hessStr()
            hessStr = sparse(tril(ones(4)));
        end
    end
end       