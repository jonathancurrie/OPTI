function test_rmathlib_formal()
% Test rmathlib against statistics toolbox
clc
clear 
tol = 1e-6;

if (~exist('tinv.m','file'))
    fprintf('The statistics toolbox is not installed, skipping test\n');
    return;
end

%% Student's t distribution
fprintf('Comparing T Distribution\n');
x = 0; q = 1.0589; p = 0.95; df = 1:10;
opticheckval.relErrorCheck(rmathlib('qt',p,df),tinv(p,df),'tinv',tol);
opticheckval.relErrorCheck(rmathlib('dt',x,df),tpdf(x,df),'tpdf',tol);
opticheckval.relErrorCheck(rmathlib('pt',q,df),tcdf(q,df),'tcdf',tol);

%% F distribution
fprintf('Comparing F Distribution\n');
x = 3; q = 3; p = 0.95; df1 = 5:15; df2 = 10:20;
opticheckval.relErrorCheck(rmathlib('qf',p,df1,df2),finv(p,df1,df2),'finv',tol);
opticheckval.relErrorCheck(rmathlib('df',x,df1,df2),fpdf(x,df1,df2),'fpdf',tol);
opticheckval.relErrorCheck(rmathlib('pf',q,df1,df2),fcdf(q,df1,df2),'fcdf',tol);

%% Normal Distribution
fprintf('Comparing Normal Distribution\n');
x = [-2 -1 0 1 2]; q = [-2,-1,0,1,2]; p = [0.1,0.25,0.5,0.75,0.9]; mu = 0; sig = 1;
opticheckval.relErrorCheck(rmathlib('qnorm',p,mu,sig),icdf('Normal',p,mu,sig),'Normal icdf',tol);
opticheckval.relErrorCheck(rmathlib('dnorm',x,mu,sig),pdf('Normal',x,mu,sig),'Normal pdf',tol);
opticheckval.relErrorCheck(rmathlib('pnorm',q,mu,sig),cdf('Normal',q,mu,sig),'Normal cdf',tol);

%% Uniform Distribution
fprintf('Comparing Uniform Distribution\n');
x = [-2 -1 0 1 2]; q = [-2,-1,0,1,2]; p = [0.1,0.25,0.5,0.75,0.9]; a = 0.1; b = 0.5;
opticheckval.relErrorCheck(rmathlib('qunif',p,a,b),icdf('Uniform',p,a,b),'Uniform icdf',tol);
opticheckval.relErrorCheck(rmathlib('dunif',x,a,b),pdf('Uniform',x,a,b),'Uniform pdf',tol);
opticheckval.relErrorCheck(rmathlib('punif',q,a,b),cdf('Uniform',q,a,b),'Uniform cdf',tol);

%% Chi2 distribution
fprintf('Comparing Chi2 Distribution\n');
x = 0; q = 1.0589; p = 0.95; df = 1:10;
opticheckval.relErrorCheck(rmathlib('qchisq',p,df),chi2inv(p,df),'chi2inv',tol);
opticheckval.relErrorCheck(rmathlib('dchisq',x,df),chi2pdf(x,df),'chi2pdf',tol);
opticheckval.relErrorCheck(rmathlib('pchisq',q,df),chi2cdf(q,df),'chi2cdf',tol);

%% Exponential distribution
fprintf('Comparing Exponential Distribution\n');
x = 0; q = 1.0589; p = 0.95; u = 1:10;
opticheckval.relErrorCheck(rmathlib('qexp',p,u),icdf('Exponential',p,u),'Exponential icdf',tol);
opticheckval.relErrorCheck(rmathlib('dexp',x,u),pdf('Exponential',x,u),'Exponential pdf',tol);
opticheckval.relErrorCheck(rmathlib('pexp',q,u),cdf('Exponential',q,u),'Exponential cdf',tol);

%% Beta Distribution
fprintf('Comparing Beta Distribution\n');
x = [-2 -1 0 1 2]; q = [-2,-1,0,1,2]; p = [0.1,0.25,0.5,0.75,0.9]; a = 0.1; b = 0.5;
opticheckval.relErrorCheck(rmathlib('qbeta',p,a,b),icdf('Beta',p,a,b),'Beta icdf',tol);
opticheckval.relErrorCheck(rmathlib('dbeta',x,a,b),pdf('Beta',x,a,b),'Beta pdf',tol);
opticheckval.relErrorCheck(rmathlib('pbeta',q,a,b),cdf('Beta',q,a,b),'Beta cdf',tol);

% %% Cauchy Distribution
% fprintf('Comparing Cauchy Distribution\n');
% x = [-2 -1 0 1 2]; q = [-2,-1,0,1,2]; p = [0.1,0.25,0.5,0.75,0.9]; a = 0.1; b = 0.5;
% opticheckval.relErrorCheck(rmathlib('qcauchy',p,a,b),icdf('Cauchy',p,a,b),'Cauchy icdf',tol);
% opticheckval.relErrorCheck(rmathlib('dcauchy',x,a,b),pdf('Cauchy',x,a,b),'Cauchy pdf',tol);
% opticheckval.relErrorCheck(rmathlib('pcauchy',q,a,b),cdf('Cauchy',q,a,b),'Cauchy cdf',tol);

%% Gamma Distribution
fprintf('Comparing Gamma Distribution\n');
x = [-2 -1 0 1 2]; q = [-2,-1,0,1,2]; p = [0.1,0.25,0.5,0.75,0.9]; a = 0.1; b = 0.5;
opticheckval.relErrorCheck(rmathlib('qgamma',p,a,b),icdf('Gamma',p,a,b),'Gamma icdf',tol);
opticheckval.relErrorCheck(rmathlib('dgamma',x,a,b),pdf('Gamma',x,a,b),'Gamma pdf',tol);
opticheckval.relErrorCheck(rmathlib('pgamma',q,a,b),cdf('Gamma',q,a,b),'Gamma cdf',tol);

%% Log Normal Distribution
fprintf('Comparing Log Normal Distribution\n');
x = [-2 -1 0 1 2]; q = [-2,-1,0,1,2]; p = [0.1,0.25,0.5,0.75,0.9]; mu = log(20000); sig = 1;
opticheckval.relErrorCheck(rmathlib('qlnorm',p,mu,sig),icdf('Lognormal',p,mu,sig),'Lognormal icdf',tol);
opticheckval.relErrorCheck(rmathlib('dlnorm',x,mu,sig),pdf('Lognormal',x,mu,sig),'Lognormal pdf',tol);
opticheckval.relErrorCheck(rmathlib('plnorm',q,mu,sig),cdf('Lognormal',q,mu,sig),'Lognormal cdf',tol);

%% Poisson distribution
fprintf('Comparing Poisson Distribution\n');
x = 0; q = 1.0589; p = 0.95; lambda = 1:10;
opticheckval.relErrorCheck(rmathlib('qpois',p,lambda),icdf('Poisson',p,lambda),'Poisson icdf',tol);
opticheckval.relErrorCheck(rmathlib('dpois',x,lambda),pdf('Poisson',x,lambda),'Poisson pdf',tol);
opticheckval.relErrorCheck(rmathlib('ppois',q,lambda),cdf('Poisson',q,lambda),'Poisson cdf',tol);

%% Weibull Distribution
fprintf('Comparing Weibull Distribution\n');
x = [-2 -1 0 1 2]; q = [-2,-1,0,1,2]; p = [0.1,0.25,0.5,0.75,0.9]; shape = 0.5; scale = 1;
opticheckval.relErrorCheck(rmathlib('qweibull',p,shape,scale),icdf('Weibull',p,scale,shape),'Weibull icdf',tol);
opticheckval.relErrorCheck(rmathlib('dweibull',x,shape,scale),pdf('Weibull',x,scale,shape),'Weibull pdf',tol);
opticheckval.relErrorCheck(rmathlib('pweibull',q,shape,scale),cdf('Weibull',q,scale,shape),'Weibull cdf',tol);
