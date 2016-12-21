function [prob,sol,fmin] = nls_prob(varargin)
%NLS_PROB  Return an OPTI NLS 
%
%   prob = nls_prob(no) return a pre-built optiprob of a saved NLS.
%
%   [prob,sol,fmin] = NLS_prob(no) returns the optimum solution and function
%   eval at the optimum
%
%   no = nls_prob() returns the number of problems available for testing.

%   (C) 2012 Jonathan Currie (IPL)

% Functions are taken from:
% More, J. J., B. S. Garbow, and K. E. Hillstrom. "Testing Unconstrained 
% Optimization Software." ACM Transactions on Mathematical Software 7, 
% no. 1 (1981): 17-41. 

%Check if just returning no problems
if(nargin < 1)
    prob = 20; sol = []; fmin = [];
    return;
else
    no = varargin{1};
end          

%Big switch yard
switch(no)
    case 1 
        fun = @(x) [10*(x(2)-x(1)^2); 1 - x(1)];
        ydata = zeros(2,1);
        x0 = [-1.2;1];
        prob = optiprob('fun',fun,'ydata',ydata,'x0',x0);            
        sol = [1,1];
        fmin = 0;
        
    case 2 
        fun = @(x) [-13 + x(1) + ((5 - x(2))*x(2) - 2)*x(2);
                    -29 + x(1) + ((x(2) + 1)*x(2) - 14)*x(2)];
        ydata = zeros(2,1);
        x0 = [0.5;-2];
        prob = optiprob('fun',fun,'ydata',ydata,'x0',x0);            
        sol = [5,4];
        fmin = 0;
        
    case 3 
        fun = @(x) [1e4*x(1)*x(2) - 1;
                    exp(-x(1)) + exp(-x(2)) - 1.0001];
        ydata = zeros(2,1);
        x0 = [0;1];
        prob = optiprob('fun',fun,'ydata',ydata,'x0',x0);            
        sol = [1.098e-5;9.106];
        fmin = 0;
        
    case 4 
        fun = @(x) [x(1) - 1e6; x(2) - 2e-6; x(1)*x(2) - 2];
        ydata = zeros(3,1);
        x0 = [1;1];
        prob = optiprob('fun',fun,'ydata',ydata,'x0',x0);            
        sol = [1e6,2e-6];
        fmin = 0;
        
    case 5 
        fun = @(x) [1.5 - x(1)*(1 - x(2)); 2.25 - x(1)*(1 - x(2)^2); 2.625 - x(1)*(1 - x(2)^3) ];
        ydata = zeros(3,1);
        x0 = [1;1];
        prob = optiprob('fun',fun,'ydata',ydata,'x0',x0);            
        sol = [3;0.5];
        fmin = 0;
        
    case 6 
        i = (1:10)';
        fun = @(x) 2 + 2*i - (exp(0.2578*i) + exp(0.2578*i));
        ydata = zeros(10,1);
        x0 = [0.3;0.4];
        prob = optiprob('fun',fun,'ydata',ydata,'x0',x0);            
        sol = [0.2578;0.2578];
        fmin = 124.36226865;
        
    case 7 
        fun = @more7;
        ydata = zeros(3,1);
        x0 = [-1;0;0];
        prob = optiprob('fun',fun,'ydata',ydata,'x0',x0);            
        sol = [1;0;0];
        fmin = 0;
        
    case 8 
        u = (1:15)';
        v = 16 - u;
        w = zeros(15,1);
        for i = 1:15
            w(i) = min(u(i),v(i));
        end
        y = [0.14;0.18;0.22;0.25;0.29;0.32;0.35;0.39;0.37;0.58;0.73;0.96;1.34;2.10;4.39];        
        fun = @(x) (x(1) + u ./ (v*x(2) + w*x(3)));
        x0 = [1;1;1];
        prob = optiprob('fun',fun,'ydata',y,'x0',x0);            
        sol = [];
        fmin = 8.21487e-3;
        
    case 9 
        i = (1:15)';
        y = [0.0009;0.0044;0.0175;0.0540;0.1295;0.2420;0.3521;0.3989;0.3521;0.2420;0.1295;0.0540;0.0175;0.0044;0.0009];
        t = (8 - i)./2;
        fun = @(x) x(1)*exp((-x(2)*(t - x(3)).^2)./2);
        x0 = [0.4;1;0];
        prob = optiprob('fun',fun,'ydata',y,'x0',x0);            
        sol = [];
        fmin = 1.12793e-8;
        
    case {10,11}
        i = (1:16)';
        y = [34780;28610;23650;19630;16370;13720;11540;9744;8261;7030;6005;5147;4427;3820;3307;2872];
        t = 45 + 5*i;
        fun = @(x) x(1)*exp(x(2)./(t + x(3)));   
        x0 = [0.02;4000;250];
        prob = optiprob('fun',fun,'ydata',y,'x0',x0);            
        sol = [];
        fmin = 87.9458;      
    
%     case 11 %not working??
%         m = 100;
%         i = [1:m]';        
%         t = i./100;
%         y = 25 + (-50 * log(t)).^(2/3);
%         fun = @(x) exp(-(abs(y*100.*i*x(2)).^x(3))./x(1)) - t; 
%         x0 = [5;2.5;0.15];
%         ydata = zeros(m,1);
%         prob = optiprob('fun',fun,'ydata',ydata,'x0',x0);            
%         sol = [50;25;1.5];
%         fmin = 0; 
        
    case 12
        m = 3;
        i = (1:m)';        
        t = 0.1*i;
        fun = @(x) exp(-t*x(1)) - exp(-t*x(2)) - x(3)*exp(-t) - exp(-10*t); 
        x0 = [0;10;20];
        ydata = zeros(m,1);
        prob = optiprob('fun',fun,'ydata',ydata,'x0',x0);            
        sol = [1;10;1];
        fmin = 0;
        
    case 13
        fun = @(x) [x(1) + 10*x(2);
                    sqrt(5)*(x(3) - x(4));
                    (x(3) - 2*x(3))^2;
                    sqrt(10)*(x(1)-x(4))^2]; 
        x0 = [3;-1;0;1];
        ydata = zeros(4,1);
        prob = optiprob('fun',fun,'ydata',ydata,'x0',x0);            
        sol = [0;0;0;0];
        fmin = 0;
        
    case 14
        fun = @(x) [10*(x(2)-x(1)^2);
                    1 - x(1);
                    sqrt(90)*(x(4)-x(3)^2);
                    1 - x(3);
                    sqrt(10)*(x(2) + x(4) - 2);
                    10^(-0.5)*(x(2)-x(4))]; 
        x0 = [3;-1;-3;-1];
        ydata = zeros(6,1);
        prob = optiprob('fun',fun,'ydata',ydata,'x0',x0);            
        sol = [1;1;1;1];
        fmin = 0;
        
    case 15
        y = [0.1957;0.1947;0.1735;0.16;0.0844;0.0627;0.0456;0.0342;0.0323;0.0235;0.0246];
        u = [4;2;1;0.5;0.25;0.167;0.125;0.1;0.0833;0.0714;0.0625];
        fun = @(x) (x(1)*(u.^2 + u*x(2))) ./ (u.^2 + u*x(3) + x(4)); 
        x0 = [0.25;0.39;0.415;0.39];
        prob = optiprob('fun',fun,'ydata',-y,'x0',x0);            
        sol = [];
        fmin = 3.07505e-4;
        
    case 16
        m = 20;
        i = (1:m)';
        t = i/5;
        fun = @(x) (x(1) + t*x(2) - exp(t)).^2 + (x(3) + x(4)*sin(t) - cos(t)).^2; 
        x0 = [25;5;-5;-1];
        ydata = zeros(m,1);
        prob = optiprob('fun',fun,'ydata',ydata,'x0',x0);            
        sol = [];
        fmin = 85822.2016;
        
    case 17
        m = 33;
        i = (1:m)';
        t = 10*(i - 1);
        y = [0.844;0.908;0.932;0.936;0.925;0.908;0.881;0.850;0.818;0.784;0.751;0.718;0.685;0.658;0.628;0.603;...
             0.58;0.558;0.538;0.522;0.506;0.490;0.478;0.467;0.457;0.448;0.438;0.431;0.424;0.42;0.414;0.411;0.406];
        fun = @(x) (x(1) + x(2)*exp(-t*x(4)) + x(3)*exp(-t*x(5))); 
        x0 = [0.5;1.5;-1;0.01;0.02];
        prob = optiprob('fun',fun,'ydata',-y,'x0',x0);            
        sol = [];
        fmin = 5.46489e-5;
        
    case 18
        m = 13;
        i = (1:m)';
        t = 0.1*i;
        y = exp(-t) - 5*exp(-10*t) + 3*exp(-4*t);
        fun = @(x) x(3)*exp(-t*x(1)) - x(4)*exp(-t*x(2)) + x(6)*exp(-t*x(5)); 
        x0 = [1;2;1;1;1;1];
        prob = optiprob('fun',fun,'ydata',y,'x0',x0);            
        sol = [];
        fmin = 0; %5.65565e-3; seems wrong?
        
    case 19
        m = 65;
        i = (1:m)';
        t = (i-1)/10;
        y = [1.366;1.191;1.112;1.013;0.991;0.885;0.831;0.847;0.786;0.725;0.746;0.679;0.608;0.655;0.616;0.606;0.602;0.626;0.651;0.724;0.649;0.649;...
            0.694;0.644;0.624;0.661;0.612;0.558;0.533;0.495;0.5;0.423;0.395;0.375;0.372;0.391;0.396;0.405;0.428;0.429;0.523;0.562;0.607;0.653;...
            0.672;0.708;0.633;0.668;0.645;0.632;0.591;0.559;0.597;0.625;0.739;0.710;0.729;0.720;0.636;0.581;0.428;0.292;0.162;0.098;0.054];
        fun = @(x) (x(1)*exp(-t*x(5)) + x(2)*exp(-(t-x(9)).^2*x(6)) + x(3)*exp(-(t-x(10)).^2*x(7)) + x(4)*exp(-(t-x(11)).^2*x(8))); 
        x0 = [1.3;0.65;0.65;0.7;0.6;3;5;7;2;4.5;5.5];
        prob = optiprob('fun',fun,'ydata',y,'x0',x0);            
        sol = [];
        fmin = 4.01377e-2;
        
    case 20
        m = 31; n = 12;
        fun = @more20; 
        x0 = zeros(n,1);
        ydata = zeros(m,1);
        prob = optiprob('fun',fun,'ydata',ydata,'x0',x0);            
        sol = [];
        fmin = 4.72238e-10;
        
    otherwise
        error('Problem not available or not implemented yet');
end


function j = more7(x)

if(x(1) > 0)
    theta = 1/(2*pi)*atan(x(2)/x(1));
else
    theta = 1/(2*pi)*atan(x(2)/x(1)) + 0.5;
end

j(1,1) = 10*(x(3) - 10*theta);
j(2,1) = 10*(x(1)^2 + x(2)^2)^0.5 - 1;
j(3,1) = x(3);


function f = more20(x)
n = length(x);
j = (2:n)';
j2 = (1:n)';
t = 0.1*(1:29)';
f = zeros(31,1);
for k = 1:29
    f(k) = sum( (j-1).*x(2:end).*t(k).^(j-2) ) - sum( x.*t(k).^(j2-1) ) - 1;
end
f(30) = x(1);
f(31) = x(2) - x(1)^2 - 1;



