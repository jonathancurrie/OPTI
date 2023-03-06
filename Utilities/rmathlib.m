% RMATHLIB  A MEX Interface to the R Standalone Math Library
%
% rmathlib uses the R Math library, An Official GNU Project (GPL v2.0)
%
%   sol = rmathlib(fcn,in1,in2,...,inn)
%
%   Input arguments:
%       fcn - R function name, see below
%       in1-inn - Function dependent input arguments, see below
%
%   Return arguments:
%       sol - solution vector
%
%
% List of R Math Functions Supported:
%
% Student's t Distribution
% - rmathlib('dt',x,df)                     - T Probability Density Function
% - rmathlib('pt',q,df)                     - T Cumulative Distribution Function
% - rmathlib('qt',p,df)                     - T Inverse CDF (Quantile Function)
% - rmathlib('rt',df,m,n)                   - T Distribution Random Number Array
%
% F Distribution
% - rmathlib('df',x,df1,df2)                - F Probability Density Function
% - rmathlib('pf',q,df1,df2)                - F Cumulative Distribution Function
% - rmathlib('qf',p,df1,df2)                - F Inverse CDF (Quantile Function)
% - rmathlib('rf',df1,df2,m,n)              - F Distribution Random Number Array
%
% Normal Distribution
% - rmathlib('dnorm',x,mu,sig)              - Normal Probability Density Function
% - rmathlib('pnorm,q,mu,sig)               - Normal Cumulative Distribution Function
% - rmathlib('qnorm,p,mu,sig)               - Normal Inverse CDF (Quantile Function)
% - rmathlib('rnorm,mu,sig,m,n)             - Normal Distribution Random Number Array
%
% Uniform Distribution
% - rmathlib('dunif,x,min,max)              - Uniform Probability Density Function
% - rmathlib('punif,q,min,max)              - Uniform Cumulative Distribution Function
% - rmathlib('qunif,p,min,max)              - Uniform Inverse CDF (Quantile Function)
% - rmathlib('runif,min,max,m,n)            - Uniform Distribution Random Number Array
%
% Chi^2 Distribution
% - rmathlib('dchisq,x,df)                  - Normal Probability Density Function
% - rmathlib('pchisq,q,df)                  - Normal Cumulative Distribution Function
% - rmathlib('qchisq,p,df)                  - Normal Inverse CDF (Quantile Function)
% - rmathlib('rchisq,df,m,n)                - Normal Distribution Random Number Array
%
% Exponential Distribution
% - rmathlib('dexp,x,rate)                  - Exponential Probability Density Function
% - rmathlib('pexp',q,rate)                 - Exponential Cumulative Distribution Function
% - rmathlib('qexp',p,rate)                 - Exponential Inverse CDF (Quantile Function)
% - rmathlib('rexp',rate,m,n)               - Exponential Distribution Random Number Array
%
% Beta Distribution
% - rmathlib('dbeta',x,shape1,shape2)       - Beta Probability Density Function
% - rmathlib('pbeta',q,shape1,shape2)       - Beta Cumulative Distribution Function
% - rmathlib('qbeta',p,shape1,shape2)       - Beta Inverse CDF (Quantile Function)
% - rmathlib('rbeta',shape1,shape2,m,n)     - Beta Distribution Random Number Array
%
% Cauchy Distribution
% - rmathlib('dcauchy',x,location,scale)    - Cauchy Probability Density Function
% - rmathlib('pcauchy',q,location,scale)    - Cauchy Cumulative Distribution Function
% - rmathlib('qcauchy',p,location,scale)    - Cauchy Inverse CDF (Quantile Function)
% - rmathlib('rcauchy',location,scale,m,n)  - Cauchy Distribution Random Number Array
%
% Gamma Distribution
% - rmathlib('dgamma',x,shape,scale)        - Gamma Probability Density Function
% - rmathlib('pgamma',q,shape,scale)        - Gamma Cumulative Distribution Function
% - rmathlib('qgamma',p,shape,scale)        - Gamma Inverse CDF (Quantile Function)
% - rmathlib('rgamma',shape,scale,m,n)      - Gamma Distribution Random Number Array
%
% Log Normal Distribution
% - rmathlib('dlnorm',x,mu_log,sig_log)     - Log Normal Probability Density Function
% - rmathlib('plnorm',q,mu_log,sig_log)     - Log Normal Cumulative Distribution Function
% - rmathlib('qlnorm',p,mu_log,sig_log)     - Log Normal Inverse CDF (Quantile Function)
% - rmathlib('rlnorm',mu_log,sig_log,m,n)   - Log Normal Distribution Random Number Array
%
% Poisson Distribution
% - rmathlib('dpois',x,lambda)              - Poisson Probability Density Function
% - rmathlib('ppois',q,lambda)              - Poisson Cumulative Distribution Function
% - rmathlib('qpois',p,lambda)              - Poisson Inverse CDF (Quantile Function)
% - rmathlib('rpois',lambda,m,n)            - Poisson Distribution Random Number Array
%
% Weibull Distribution
% - rmathlib('dweibull',x,shape,scale)      - Weibull Probability Density Function
% - rmathlib('pweibull',q,shape,scale)      - Weibull Cumulative Distribution Function
% - rmathlib('qweibull',p,shape,scale)      - Weibull Inverse CDF (Quantile Function)
% - rmathlib('rweibull',shape,scale,m,n)    - Weibull Distribution Random Number Array        

%   Copyright (C) 2017 Jonathan Currie (Control Engineering)