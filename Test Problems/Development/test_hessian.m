%% Hessian Testing
% Tests to confirm Hessian is handled correctly by OPTI

% NOTE in all cases the Hessian is the Hessian of the Lagrangian:

% sigma*H + lambda(1)*G1 + ... + lambda(M)*GM

% This means the Hessian contains the second derivatives of both the
% objective and the constraints!

%% HS71 Normal Problem IPOPT
clc
clear
f = @(x) x(1)*x(4)*sum(x(1:3)) + x(3);
g = @(x) [ x(1)*x(4) + x(4)*sum(x(1:3))
                x(1)*x(4)
                x(1)*x(4) + 1
                x(1)*sum(x(1:3)) ]';
c = @(x) [ prod(x); sum(x.^2) ];      
j = @(x) [ prod(x)./x'; 2*x' ];

cl = [25;40];
cu = [Inf;40];

H = @(x,sigma,lambda) sigma*[ 2*x(4)             0      0   0;
                              x(4)               0      0   0;
                              x(4)               0      0   0;
                              2*x(1)+x(2)+x(3)  x(1)  x(1)  0 ] + ...
                   lambda(1)*[   0          0         0         0;
                              x(3)*x(4)     0         0         0;
                              x(2)*x(4) x(1)*x(4)     0         0;
                              x(2)*x(3) x(1)*x(3) x(1)*x(2)     0  ] + ...
         		   lambda(2)*diag([2 2 2 2]);

lb = ones(4,1);
ub = 5*ones(4,1);        

x0 = [1 5 5 1]';
opts = optiset('display','iter','derivCheck','on');               

Opt = opti('fun',f,'grad',g,'hess',H,'nl',c,cl,cu,'nljac',j,'bounds',lb,ub,'options',opts)  

[x,fval,exitflag,info]= solve(Opt,x0)

%% HS71 Sense Test
clc
fm = @(x) -f(x);
gm = @(x) -g(x);
Hm = @(x,sigma,lambda) H(x,-sigma,lambda);

Opt = opti('fun',fm,'grad',gm,'Hess',Hm,'nl',c,cl,cu,'nljac',j,'bounds',lb,ub,'options',opts,'sense',-1)  

[x,fval,exitflag,info]= solve(Opt,x0)

%% HS71 MATLAB
clc
nlrhs = [25;40];
nle = [1;0];
opts = optiset('display','iter','solver','matlab');  
Opt = opti('fun',f,'grad',g,'Hess',H,'nlmix',c,nlrhs,nle,'nljac',j,'bounds',lb,ub,'options',opts)  

[x,fval,exitflag,info]= solve(Opt,x0)


%% HS71 Hessian Manual Testing
clc
x = [1.99;4.01;4.01;1.99]
lambda = [0.185942663470295
          -0.946138746645349];
%Our function
H1 = H(x,1,lambda)      

%Numdiff (uses John's code)
if(exist('derivest','file'))
    h1 = hessian(f,x)
    h2 = hessian(@(x) prod(x),x)
    h3 = hessian(@(x) sum(x.^2),x)
    H2 = h1+lambda(1)*h2+lambda(2)*h3
else
    optiwarn('OPTI:HessianTest','Derivest suite not found, no manual Hessian Tests can be checked');
end


%% HS71 Manual fmincon
clc
if (exist('fmincon.m','file'))
    mopts = optimset('display','iter','Algorithm','interior-point','GradObj','on','GradConstr','on','Hessian','user-supplied','HessFcn',@hs71H);

    [x,fv,e,i] = fmincon(@hs71F,x0,[],[],[],[],lb,ub,@hs71C,mopts)
end


%% HS71 From AMPL Comparison
clc
x = [1.99;4.01;4.01;1.99];
sigma = 0.5;
lambda = [0.185942663470295
          -0.946138746645349];
      
p = which('hs71.nl');

%Read in Problem Data
asl('open',p);

fun     = @(x) asl('fun',x);
grad    = @(x) asl('grad',x);
con     = @(x) asl('con',x);
jac     = @(x) asl('jac',x,1); %get sparse jac
jacstr  = @() asl('jacstr');
hess    = @(x,sigma,lambda) asl('hess',x,sigma,lambda);
hstr    = @() asl('hesstr');

f1 = fun(x)
f2 = f(x)

g1 = grad(x)
g2 = g(x)

c1 = con(x)
c2 = c(x)

j1 = full(jac(x))
j2 = j(x)

h1 = full(tril(hess(x,sigma,lambda)))
h2 = H(x,sigma,lambda)


asl('close');

%% HS71 From AMPL Comparison, No Grad Call
clc
x = [1.99;4.01;4.01;1.99];
sigma = 0.5;
lambda = [0.185942663470295
          -0.946138746645349];
      
p = which('hs71.nl');

%Read in Problem Data
asl('open',p);

fun     = @(x) asl('fun',x);
grad    = @(x) asl('grad',x);
con     = @(x) asl('con',x);
hess    = @(x,sigma,lambda) asl('hess',x,sigma,lambda);

% f1 = fun(x)
% f2 = f(x)
% 
% g1 = grad(x)
% g2 = g(x)

% c1 = con(x)
% c2 = c(x)

h1 = full(tril(hess(x,sigma,lambda)))
h2 = H(x,sigma,lambda)


asl('close');

%% Complete amplRead HS71 test IPOPT
clc
clear
prob = amplRead('hs71');
opts = optiset('display','iter','derivCheck','on');

O = opti(prob,opts)

[xx,ff,ee,ii] = solve(O)

%% Complete amplRead HS71 test MATLAB
clc
clear
prob = amplRead('hs71');
opts = optiset('display','iter','solver','matlab');

O = opti(prob,opts)

[xx,ff,ee,ii] = solve(O)