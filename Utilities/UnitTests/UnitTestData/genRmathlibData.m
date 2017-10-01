%% Generate unit test data for rmathlib from MATLAB's statistics toolbox
% J.Currie October 2017
clc
clear

%% Student's t Distribution
t.x       = 0;
t.q       = 1.0589;
t.p       = 0.95;
t.df      = 1:10;
t.qt      = tinv(t.p,t.df);
t.dt      = tpdf(t.x,t.df);
t.pt      = tcdf(t.q,t.df);

%% F Distribution
f.x     = 3; 
f.q     = 3; 
f.p     = 0.95; 
f.df1   = 5:15; 
f.df2   = 10:20;
f.qf    = finv(f.p,f.df1,f.df2);
f.df    = fpdf(f.x,f.df1,f.df2);
f.pf    = fcdf(f.q,f.df1,f.df2);

%% Normal Distribution
norm.x      = [-2 -1 0 1 2]; 
norm.q      = [-2,-1,0,1,2]; 
norm.p      = [0.1,0.25,0.5,0.75,0.9]; 
norm.mu     = 0; 
norm.sig    = 1;
norm.qnorm  = icdf('Normal',norm.p,norm.mu,norm.sig);
norm.dnorm  = pdf('Normal',norm.x,norm.mu,norm.sig);
norm.pnorm  = cdf('Normal',norm.q,norm.mu,norm.sig);

%% Uniform Distribution
unif.x      = [-2 -1 0 1 2]; 
unif.q      = [-2,-1,0,1,2]; 
unif.p      = [0.1,0.25,0.5,0.75,0.9]; 
unif.a      = 0.1; 
unif.b      = 0.5;
unif.qunif  = icdf('Uniform',unif.p,unif.a,unif.b);
unif.dunif	= pdf('Uniform',unif.x,unif.a,unif.b);
unif.punif	= cdf('Uniform',unif.q,unif.a,unif.b);

%% Chi2 distribution
chi2.x      = 0; 
chi2.q      = 1.0589; 
chi2.p      = 0.95; 
chi2.df     = 1:10;
chi2.qchisq	= chi2inv(chi2.p,chi2.df);
chi2.dchisq	= chi2pdf(chi2.x,chi2.df);
chi2.pchisq	= chi2cdf(chi2.q,chi2.df);

%% Exponential distribution
exp.x       = 0; 
exp.q       = 1.0589; 
exp.p       = 0.95; 
exp.u       = 1:10;
exp.qexp	= icdf('Exponential',exp.p,exp.u);
exp.dexp	= pdf('Exponential',exp.x,exp.u);
exp.pexp	= cdf('Exponential',exp.q,exp.u);

%% Beta Distribution
beta.x      = [-2 -1 0 1 2]; 
beta.q      = [-2,-1,0,1,2]; 
beta.p      = [0.1,0.25,0.5,0.75,0.9]; 
beta.a      = 0.1; 
beta.b      = 0.5;
beta.qbeta	= icdf('Beta',beta.p,beta.a,beta.b);
beta.dbeta	= pdf('Beta',beta.x,beta.a,beta.b);
beta.pbeta	= cdf('Beta',beta.q,beta.a,beta.b);

%% Cauchy Distribution
% cauchy.x    = [-2 -1 0 1 2]; 
% cauchy.q    = [-2,-1,0,1,2]; 
% cauchy.p    = [0.1,0.25,0.5,0.75,0.9]; 
% cauchy.a    = 0.1; 
% cauchy.b    = 0.5;
% cauchy.qcauchy = icdf('Cauchy',cauchy.p,cauchy.a,cauchy.b);
% cauchy.dcauchy = pdf('Cauchy',cauchy.x,cauchy.a,cauchy.b);
% cauchy.pcauchy = cdf('Cauchy',cauchy.q,cauchy.a,cauchy.b);

%% Gamma Distribution
gam.x       = [-2 -1 0 1 2]; 
gam.q       = [-2,-1,0,1,2]; 
gam.p       = [0.1,0.25,0.5,0.75,0.9]; 
gam.a       = 0.1; 
gam.b       = 0.5;
gam.qgamma	= icdf('Gamma',gam.p,gam.a,gam.b);
gam.dgamma	= pdf('Gamma',gam.x,gam.a,gam.b);
gam.pgamma	= cdf('Gamma',gam.q,gam.a,gam.b);

%% Log Normal Distribution
logn.x      = [-2 -1 0 1 2]; 
logn.q      = [-2,-1,0,1,2]; 
logn.p      = [0.1,0.25,0.5,0.75,0.9]; 
logn.mu     = log(20000); 
logn.sig    = 1;
logn.qlnorm	= icdf('Lognormal',logn.p,logn.mu,logn.sig);
logn.dlnorm = pdf('Lognormal',logn.x,logn.mu,logn.sig);
logn.plnorm = cdf('Lognormal',logn.q,logn.mu,logn.sig);

%% Poisson distribution
pois.x      = 0; 
pois.q      = 1.0589; 
pois.p      = 0.95; 
pois.lambda = 1:10;
pois.qpois  = icdf('Poisson',pois.p,pois.lambda);
pois.dpois	= pdf('Poisson',pois.x,pois.lambda);
pois.ppois	= cdf('Poisson',pois.q,pois.lambda);

%% Weibull Distribution
weib.x      = [-2 -1 0 1 2]; 
weib.q      = [-2,-1,0,1,2]; 
weib.p      = [0.1,0.25,0.5,0.75,0.9]; 
weib.shape  = 0.5; 
weib.scale  = 1;
weib.qweibull = icdf('Weibull',weib.p,weib.scale,weib.shape);
weib.dweibull = pdf('Weibull',weib.x,weib.scale,weib.shape);
weib.pweibull = cdf('Weibull',weib.q,weib.scale,weib.shape);

%% Matrix Tests
mat.tx    = 0;
mat.tq    = 1.0589;
mat.tdf   = magic(3);
mat.dt    = tpdf(mat.tx,mat.tdf);
mat.pt    = tcdf(mat.tq,mat.tdf);

mat.nx     = 1; 
mat.nq     = 1; 
mat.mu     = 0; 
mat.sig    = magic(3);
mat.dnorm  = pdf('Normal',mat.nx,mat.mu,mat.sig);
mat.pnorm  = cdf('Normal',mat.nq,mat.mu,mat.sig);

%% Upper Tests
upper.tq    = 0.1;
upper.tdf   = 0.1;
upper.pt    = tcdf(upper.tq, upper.tdf, 'upper');

upper.nq    = 0.1;
upper.mu    = 0;
upper.sig   = 1;
upper.pnorm = cdf('Normal',upper.nq, upper.mu, upper.sig, 'upper');


%% Save to mat
save Utilities/UnitTests/UnitTestData/rmathlibData