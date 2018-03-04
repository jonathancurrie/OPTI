%% A Unit Test Class for rmathlib
classdef rmathlib_tests < matlab.unittest.TestCase

    properties
        absTol = 1e-10;
        data   = load('Utilities/UnitTests/UnitTestData/rmathlibData.mat');
    end
    
    % Unit Tests
    methods (Test)
        
        %-- Input Args --%
        function inputArgs(testCase)
            testCase.verifyError(@() rmathlib(1), 'OPTIMex:InputError');   % not enough args
            testCase.verifyError(@() rmathlib(1,2), 'OPTIMex:InputError');   % not enough args
            testCase.verifyError(@() rmathlib(1,2,3), 'OPTIMex:InputError');   % not string
            testCase.verifyError(@() rmathlib('a','a',3), 'OPTIMex:InputError');   % not double
            testCase.verifyError(@() rmathlib('a',2,'a'), 'OPTIMex:InputError');   % not double
            testCase.verifyError(@() rmathlib('a',2,3), 'OPTIMex:InputError');   % unknown function
            % D1
            testCase.verifyError(@() rmathlib('rt',1,2), 'OPTIMex:InputError');     % not enough args
            testCase.verifyError(@() rmathlib('rt','a',2,3), 'OPTIMex:InputError'); % not double
            testCase.verifyError(@() rmathlib('rt',1,0,3), 'OPTIMex:InputError');   % not in bounds
            testCase.verifyError(@() rmathlib('rt',1,2,0), 'OPTIMex:InputError');   % not in bounds
            testCase.verifyError(@() rmathlib('rt',[1;1],2,3), 'OPTIMex:InputError'); % not scalar
            testCase.verifyError(@() rmathlib('rt',1,[2;2],3), 'OPTIMex:InputError'); % not scalar
            testCase.verifyError(@() rmathlib('rt',1,2,[3;3]), 'OPTIMex:InputError'); % not scalar
            % D2
            testCase.verifyError(@() rmathlib('rnorm',1,2,3), 'OPTIMex:InputError');   % not enough args
            testCase.verifyError(@() rmathlib('rnorm','a',2,3,4), 'OPTIMex:InputError'); % not double
            testCase.verifyError(@() rmathlib('rnorm',1,'a',3,4), 'OPTIMex:InputError'); % not double
            testCase.verifyError(@() rmathlib('rnorm',1,2,0,4), 'OPTIMex:InputError');   % not in bounds
            testCase.verifyError(@() rmathlib('rnorm',1,2,3,0), 'OPTIMex:InputError');   % not in bounds 
            testCase.verifyError(@() rmathlib('rnorm',[1;1],2,3,4), 'OPTIMex:InputError'); % not scalar
            testCase.verifyError(@() rmathlib('rnorm',1,[2;2],3,4), 'OPTIMex:InputError'); % not scalar
            testCase.verifyError(@() rmathlib('rnorm',1,2,[3;3],4), 'OPTIMex:InputError'); % not scalar
            testCase.verifyError(@() rmathlib('rnorm',1,2,3,[4;4]), 'OPTIMex:InputError'); % not scalar
            % D2I1
            testCase.verifyError(@() rmathlib('dt','a',2), 'OPTIMex:InputError'); % not double
            testCase.verifyError(@() rmathlib('dt',1,'a'), 'OPTIMex:InputError'); % not double
            % D2I2
            testCase.verifyError(@() rmathlib('pt','a',2), 'OPTIMex:InputError'); % not double
            testCase.verifyError(@() rmathlib('pt',1,'a'), 'OPTIMex:InputError'); % not double
            % D3I1
            testCase.verifyError(@() rmathlib('df','a',2,3), 'OPTIMex:InputError'); % not double
            testCase.verifyError(@() rmathlib('df',1,'a',3), 'OPTIMex:InputError'); % not double
            testCase.verifyError(@() rmathlib('df',1,2,'a'), 'OPTIMex:InputError'); % not double
            % D3I2
            testCase.verifyError(@() rmathlib('pf','a',2,3), 'OPTIMex:InputError'); % not double
            testCase.verifyError(@() rmathlib('pf',1,'a',3), 'OPTIMex:InputError'); % not double
            testCase.verifyError(@() rmathlib('pf',1,2,'a'), 'OPTIMex:InputError'); % not double
            % Scalar expansion checked separately below  
            
            % Check case-insensitive
            testCase.verifyEqual(rmathlib('dt',0,1), rmathlib('dT',0,1));
            testCase.verifyEqual(rmathlib('dt',0,1), rmathlib('Dt',0,1));
            testCase.verifyEqual(rmathlib('dt',0,1), rmathlib('DT',0,1));
        end
        
        %-- Scalar Expansion --%
        function validateScalarExp(testCase)
            % D2I1 - dt
            testCase.verifyEqual(1, numel(rmathlib('dt',0,1)));
            testCase.verifyEqual(2, numel(rmathlib('dt',[0;1],1)));  
            testCase.verifyEqual(3, numel(rmathlib('dt',1,[0;1;2]))); 
            testCase.verifyError(@() rmathlib('dt',[0;1],[0;1;2]), 'OPTIMex:InputError');
            testCase.verifyError(@() rmathlib('dt',[0;1;2],[0;1]), 'OPTIMex:InputError');       
            % D2I2 - pt
            testCase.verifyEqual(1, numel(rmathlib('pt',0,1)));
            testCase.verifyEqual(2, numel(rmathlib('pt',[0;1],1)));  
            testCase.verifyEqual(3, numel(rmathlib('pt',1,[0;1;2]))); 
            testCase.verifyError(@() rmathlib('pt',[0;1],[0;1;2]), 'OPTIMex:InputError');
            testCase.verifyError(@() rmathlib('pt',[0;1;2],[0;1]), 'OPTIMex:InputError'); 
            % D3I1 - dnorm
            testCase.verifyEqual(1, numel(rmathlib('dnorm',0,1,2)));
            testCase.verifyEqual(2, numel(rmathlib('dnorm',[0;1],1,1)));  
            testCase.verifyEqual(3, numel(rmathlib('dnorm',1,[0;1;2],1))); 
            testCase.verifyEqual(4, numel(rmathlib('dnorm',1,1,[0;1;2;3]))); 
            testCase.verifyEqual(2, numel(rmathlib('dnorm',[0;1],[0;1],1)));
            testCase.verifyEqual(2, numel(rmathlib('dnorm',[0;1],1,[0;1])));
            testCase.verifyEqual(2, numel(rmathlib('dnorm',1,[0;1],[0;1])));
            testCase.verifyError(@() rmathlib('dnorm',1,[0;1],[0;1;2]), 'OPTIMex:InputError'); 
            testCase.verifyError(@() rmathlib('dnorm',1,[0;1;2],[0;1]), 'OPTIMex:InputError');             
            testCase.verifyError(@() rmathlib('dnorm',[0;1],[0;1;2],1), 'OPTIMex:InputError'); 
            testCase.verifyError(@() rmathlib('dnorm',[0;1;2],[0;1],1), 'OPTIMex:InputError');             
            testCase.verifyError(@() rmathlib('dnorm',[0;1],1,[0;1;2]), 'OPTIMex:InputError'); 
            testCase.verifyError(@() rmathlib('dnorm',[0;1;2],1,[0;1]), 'OPTIMex:InputError');             
            testCase.verifyError(@() rmathlib('dnorm',[0;1],[0;1],[0;1;2]), 'OPTIMex:InputError'); 
            testCase.verifyError(@() rmathlib('dnorm',[0;1;2],[0;1],[0;1]), 'OPTIMex:InputError'); 
            testCase.verifyError(@() rmathlib('dnorm',[0;1],[0;1;2],[0;1]), 'OPTIMex:InputError');          
            % D3I2 - pnorm
            testCase.verifyEqual(1, numel(rmathlib('pnorm',0,1,2)));
            testCase.verifyEqual(2, numel(rmathlib('pnorm',[0;1],1,1)));  
            testCase.verifyEqual(3, numel(rmathlib('pnorm',1,[0;1;2],1))); 
            testCase.verifyEqual(4, numel(rmathlib('pnorm',1,1,[0;1;2;3]))); 
            testCase.verifyEqual(2, numel(rmathlib('pnorm',[0;1],[0;1],1)));
            testCase.verifyEqual(2, numel(rmathlib('pnorm',[0;1],1,[0;1])));
            testCase.verifyEqual(2, numel(rmathlib('pnorm',1,[0;1],[0;1])));
            testCase.verifyError(@() rmathlib('pnorm',1,[0;1],[0;1;2]), 'OPTIMex:InputError'); 
            testCase.verifyError(@() rmathlib('pnorm',1,[0;1;2],[0;1]), 'OPTIMex:InputError');             
            testCase.verifyError(@() rmathlib('pnorm',[0;1],[0;1;2],1), 'OPTIMex:InputError'); 
            testCase.verifyError(@() rmathlib('pnorm',[0;1;2],[0;1],1), 'OPTIMex:InputError');             
            testCase.verifyError(@() rmathlib('pnorm',[0;1],1,[0;1;2]), 'OPTIMex:InputError'); 
            testCase.verifyError(@() rmathlib('pnorm',[0;1;2],1,[0;1]), 'OPTIMex:InputError');             
            testCase.verifyError(@() rmathlib('pnorm',[0;1],[0;1],[0;1;2]), 'OPTIMex:InputError'); 
            testCase.verifyError(@() rmathlib('pnorm',[0;1;2],[0;1],[0;1]), 'OPTIMex:InputError'); 
            testCase.verifyError(@() rmathlib('pnorm',[0;1],[0;1;2],[0;1]), 'OPTIMex:InputError');                         
        end
        
        %-- Matrix Test --%
        function validateMatrix(testCase)
            mat = testCase.data.mat;
            % D2I1 - dt
            testCase.verifyEqual(mat.dt, rmathlib('dt',mat.tx,mat.tdf), 'AbsTol', testCase.absTol);    
            % D2I2 - pt
            testCase.verifyEqual(mat.pt, rmathlib('pt',mat.tq,mat.tdf), 'AbsTol', testCase.absTol);
            % D3I1 - dnorm
            testCase.verifyEqual(mat.dnorm, rmathlib('dnorm',mat.nx,mat.mu,mat.sig), 'AbsTol', testCase.absTol);     
            % D3I2 - pnorm
            testCase.verifyEqual(mat.pnorm, rmathlib('pnorm',mat.nq,mat.mu,mat.sig), 'AbsTol', testCase.absTol);
        end
        
        %-- Integer Args test --%
        function validateIntegerArgs(testCase)
            upper = testCase.data.upper;
            % D2I1 - dt
            dt = rmathlib('dt',0,1);
            testCase.verifyEqual(log(dt), rmathlib('dt',0,1,true), 'AbsTol', testCase.absTol);    
            % D2I2 - pt
            ptDef           = rmathlib('pt',0.1,0.1);
            ptLowerNoLog    = rmathlib('pt',0.1,0.1,true,false);
            ptUpperNoLog    = rmathlib('pt',0.1,0.1,false,false);
            ptLowerLog      = rmathlib('pt',0.1,0.1,true,true);
            ptUpperLog      = rmathlib('pt',0.1,0.1,false,true);
            testCase.verifyEqual(ptLowerNoLog, ptDef, 'AbsTol', testCase.absTol); 
            testCase.verifyEqual(log(ptLowerNoLog), ptLowerLog, 'AbsTol', testCase.absTol); 
            testCase.verifyEqual(log(ptUpperNoLog), ptUpperLog, 'AbsTol', testCase.absTol); 
            testCase.verifyEqual(upper.pt, rmathlib('pt',upper.tq, upper.tdf, false), 'AbsTol', testCase.absTol); 
            % D3I1 - dnorm
            dnorm = rmathlib('dnorm',0,0,1);
            testCase.verifyEqual(log(dnorm), rmathlib('dnorm',0,0,1,true), 'AbsTol', testCase.absTol);     
            % D3I2 - pnorm
            pnormDef           = rmathlib('pnorm',0.1,0,1);
            pnormLowerNoLog    = rmathlib('pnorm',0.1,0,1,true,false);
            pnormUpperNoLog    = rmathlib('pnorm',0.1,0,1,false,false);
            pnormLowerLog      = rmathlib('pnorm',0.1,0,1,true,true);
            pnormUpperLog      = rmathlib('pnorm',0.1,0,1,false,true);
            testCase.verifyEqual(pnormLowerNoLog, pnormDef, 'AbsTol', testCase.absTol); 
            testCase.verifyEqual(log(pnormLowerNoLog), pnormLowerLog, 'AbsTol', testCase.absTol); 
            testCase.verifyEqual(log(pnormUpperNoLog), pnormUpperLog, 'AbsTol', testCase.absTol); 
            testCase.verifyEqual(upper.pnorm, rmathlib('pnorm',upper.nq, upper.mu, upper.sig, false), 'AbsTol', testCase.absTol); 
        end
        
        
        %-- Student's t distribution --%
        function validateStudent(testCase)
            t = testCase.data.t;
            testCase.verifyEqual(t.qt, rmathlib('qt',t.p,t.df), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(t.dt, rmathlib('dt',t.x,t.df), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(t.pt, rmathlib('pt',t.q,t.df), 'AbsTol', testCase.absTol);
        end
        
        %-- F distribution --%
        function validateF(testCase)
            f = testCase.data.f;
            testCase.verifyEqual(f.qf, rmathlib('qf',f.p,f.df1,f.df2), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(f.df, rmathlib('df',f.x,f.df1,f.df2), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(f.pf, rmathlib('pf',f.q,f.df1,f.df2), 'AbsTol', testCase.absTol);
        end
        
        %-- Normal distribution --%
        function validateNormal(testCase)
            norm = testCase.data.norm;
            testCase.verifyEqual(norm.qnorm, rmathlib('qnorm',norm.p,norm.mu,norm.sig), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(norm.dnorm, rmathlib('dnorm',norm.x,norm.mu,norm.sig), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(norm.pnorm, rmathlib('pnorm',norm.q,norm.mu,norm.sig), 'AbsTol', testCase.absTol);
        end
        
        %-- Uniform distribution --%
        function validateUniform(testCase)
            unif = testCase.data.unif;
            testCase.verifyEqual(unif.qunif, rmathlib('qunif',unif.p,unif.a,unif.b), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(unif.dunif, rmathlib('dunif',unif.x,unif.a,unif.b), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(unif.punif, rmathlib('punif',unif.q,unif.a,unif.b), 'AbsTol', testCase.absTol);
        end
        
        %-- Chi2 distribution --%
        function validateChi2(testCase)
            chi2 = testCase.data.chi2;
            testCase.verifyEqual(chi2.qchisq, rmathlib('qchisq',chi2.p,chi2.df), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(chi2.dchisq, rmathlib('dchisq',chi2.x,chi2.df), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(chi2.pchisq, rmathlib('pchisq',chi2.q,chi2.df), 'AbsTol', testCase.absTol);
        end
        
        %-- Exponential distribution --%
        function validateExponential(testCase)
            exp = testCase.data.exp;
            testCase.verifyEqual(exp.qexp, rmathlib('qexp',exp.p,exp.u), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(exp.dexp, rmathlib('dexp',exp.x,exp.u), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(exp.pexp, rmathlib('pexp',exp.q,exp.u), 'AbsTol', testCase.absTol);
        end
        
        %-- Beta distribution --%
        function validateBeta(testCase)
            beta = testCase.data.beta;
            testCase.verifyEqual(beta.qbeta, rmathlib('qbeta',beta.p,beta.a,beta.b), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(beta.dbeta, rmathlib('dbeta',beta.x,beta.a,beta.b), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(beta.pbeta, rmathlib('pbeta',beta.q,beta.a,beta.b), 'AbsTol', testCase.absTol);
        end
        
        %-- Gamma distribution --%
        function validateGamma(testCase)
            gam = testCase.data.gam;
            testCase.verifyEqual(gam.qgamma, rmathlib('qgamma',gam.p,gam.a,gam.b), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(gam.dgamma, rmathlib('dgamma',gam.x,gam.a,gam.b), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(gam.pgamma, rmathlib('pgamma',gam.q,gam.a,gam.b), 'AbsTol', testCase.absTol);
        end
        
        %-- Log-Normal distribution --%
        function validateLogNormal(testCase)
            logn = testCase.data.logn;
            testCase.verifyEqual(logn.qlnorm, rmathlib('qlnorm',logn.p,logn.mu,logn.sig), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(logn.dlnorm, rmathlib('dlnorm',logn.x,logn.mu,logn.sig), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(logn.plnorm, rmathlib('plnorm',logn.q,logn.mu,logn.sig), 'AbsTol', testCase.absTol);
        end
        
        %-- Poisson distribution --%
        function validatePoisson(testCase)
            pois = testCase.data.pois;
            testCase.verifyEqual(pois.qpois, rmathlib('qpois',pois.p,pois.lambda), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(pois.dpois, rmathlib('dpois',pois.x,pois.lambda), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(pois.ppois, rmathlib('ppois',pois.q,pois.lambda), 'AbsTol', testCase.absTol);
        end
        
        %-- Weibull distribution --%
        function validateWeibull(testCase)
            weib = testCase.data.weib;
            testCase.verifyEqual(weib.qweibull, rmathlib('qweibull',weib.p,weib.shape,weib.scale), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(weib.dweibull, rmathlib('dweibull',weib.x,weib.shape,weib.scale), 'AbsTol', testCase.absTol);
            testCase.verifyEqual(weib.pweibull, rmathlib('pweibull',weib.q,weib.shape,weib.scale), 'AbsTol', testCase.absTol);
        end
    end
    
end

