function [prob,sol,fmin] = lp_prob(varargin)
%LP_PROB  Return an OPTI LP 
%
%   prob = lp_prob(no) return a pre-built optiprob of a saved LP.
%
%   [prob,sol,fmin] = lp_prob(no) returns the optimum solution and function
%   eval at the optimum
%
%   no = lp_prob() returns the number of problems available for testing.

%   (C) 2011 Jonathan Currie (IPL)

%Check if just returning no problems
if(nargin < 1)
    prob = 11; sol = []; fmin = [];
    return;
else
    no = varargin{1};
end          

%Big switch yard
switch(no)
    case 1 
        f = -[-1, 2]';
        A = [2, 1;-4, 4];
        b = [5, 5]';
        e = -[1, 1];    
        prob = optiprob('f',f,'mix',A,b,e,'name','TestLP1');            
        sol = [1.25;2.5];
        fmin = -3.75;
        
    case 2 
        f = -[50, 100];
        A = [10, 5;4, 10; 1, 1.5];
        b = [2500, 2000, 450]';
        e = [-1, -1, -1];    
        prob = optiprob('f',f,'mix',A,b,e,'name','TestLP2');            
        sol = [187.5;125];
        fmin = -21875;
        
    case 3 
        f = [40, 36];
        A = [5, 3];
        b = 45;
        e = 1;
        ub = [8, 10];    
        prob = optiprob('f',f,'mix',A,b,e,'ub',ub,'name','TestLP3');            
        sol = [8;5/3];
        fmin = 380;
        
    case 4 
        f = [3, -7, -12];
        A = [-3, 6, 8;6, -3, 7;-6, 3, 3];
        b = [12, 8, 5]';
        e = [-1, -1, -1];   
        prob = optiprob('f',f,'mix',A,b,e,'name','TestLP4');            
        sol = [-0.0666666666;0.23333333333;1.3];
        fmin = -17.43333333333333;
        
    case 5 
        f = [2, 3, 7, 7];
        A = [1, 1, -2, -5;-1, 2, 1, 4];
        b = [2, -3]';
        e = [1, 1];
        lb = zeros(4,1); 
        ub = [30 100 20 1]';   
        prob = optiprob('f',f,'mix',A,b,e,'bounds',lb,ub,'name','TestLP5');            
        sol = [2;0;0;0];
        fmin = 4;
        
    case 6 
        f = [1, 2, 3, 7, 8, 8];
        A = [5, -3, 2, -3, -1, 2; -1, 0, 2, 1, 3, -3;1, 2, -1, 0, 5, -1];
        b = [-5, -1, 3]';
        e = [1, 1, 1];
        lb = zeros(6,1);
        ub = 10*ones(6,1);   
        prob = optiprob('f',f,'mix',A,b,e,'bounds',lb,ub,'name','TestLP6');            
        sol = [0;1.5;0;0;0;0];
        fmin = 3;
        
    case 7 
        n = 40;
        t = (0:n-1)';
        y = 3.5 -.2*t;
        b = y + 0.5*ones(size(y));
        m = [ones(n,1),t(:)];
        A = [m,-m,eye(n)];
        f = [sum(m),sum(-m),2*ones(1,n)];
        e = ones(n,1);
        lb = zeros(n+4,1);
        ub = [10, 10, 10, 10, 5*ones(1,n)];  
        prob = optiprob('f',f,'mix',A,b,e,'bounds',lb,ub,'name','TestLP7');            
        sol = [4;0;0;0.2;zeros(40,1)];
        fmin = 4.0000000000000044;
        
    case 8 
        f = -[8, 15];
        A = [10, 21;2, 1];
        b = [156, 22]';
        e = [-1, -1]; 
        prob = optiprob('f',f,'mix',A,b,e,'name','TestLP8');            
        sol = [9.5625;2.875];
        fmin = -119.625;
        
    case 9 
        f = -[3, 13];
        A = [2, 9;11, -8];
        b = [40, 82]';
        e = [-1, -1]; 
        prob = optiprob('f',f,'mix',A,b,e,'name','TestLP9');            
        sol = [9.2;2.4];
        fmin = -58.8;
        
    case 10 
        f = -[592, 381, 273, 55, 48, 37, 23];
        A = [3534, 2356, 1767, 589, 528, 451, 304];
        b = 119567;
        e = -1;
        lb = zeros(7,1); 
        ub = [100 50 33 20 77 44 20]';
        prob = optiprob('f',f,'mix',A,b,e,'bounds',lb,ub,'name','TestLP10');            
        sol = [33.8333333333333;zeros(6,1)];
        fmin = -20029.33333333333;        
        
    case 11
        f = [-1 -1 -3 -2 -2]';
        A = [-1 -1 1 1 0;
             1 0 1 -3 0];
        b = [30;30];
        ub = [40;1;inf;inf;1];
        prob = optiprob('f',f,'ineq',A,b,'ub',ub,'name','TestLP11');            
        sol = [40;1;50.75;20.25;1];
        fmin = -235.7500; 
        
    otherwise
        error('Problem not available or not implemented yet');
end