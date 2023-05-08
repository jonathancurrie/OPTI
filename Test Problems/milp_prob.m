function [prob,sol,fmin] = milp_prob(varargin)
%MILP_PROB  Return an OPTI MILP 
%
%   prob = milp_prob(no) return a pre-built optiprob of a saved MILP.
%
%   [prob,sol,fmin] = milp_prob(no) returns the optimum solution and 
%   function eval at the optimum
%
%   no = milp_prob() returns the number of problems available for testing.

%   (C) 2011 Jonathan Currie (Control Engineering)

%Check if just returning no problems
if(nargin < 1)
    prob = 10; sol = []; fmin = [];
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
        xint = 'II';
        prob = optiprob('f',f,'mix',A,b,e,'int',xint);            
        sol = [1;2];
        fmin = -3;
        
    case 2 
        f = -[50, 100];
        A = [10, 5;4, 10; 1, 1.5];
        b = [2500, 2000, 450]';
        e = [-1, -1, -1];  
        xint = 'II';
        prob = optiprob('f',f,'mix',A,b,e,'int',xint);            
        sol = [187;125];
        fmin = -21850;
        
    case 3 
        f = [40, 36];
        A = [5, 3];
        b = 45;
        e = 1;
        ub = [8, 10]; 
        xint = 'II';
        prob = optiprob('f',f,'mix',A,b,e,'ub',ub,'int',xint);            
        sol = [8;2];
        fmin = 392;
        
    case 4 
        f = [3, -7, -12];
        A = [-3, 6, 8;6, -3, 7;-6, 3, 3];
        b = [12, 8, 5]';
        e = [-1, -1, -1];   
        xint = 'III';
        prob = optiprob('f',f,'mix',A,b,e,'int',xint);            
        sol = [2;3;0];
        fmin = -15;
        
    case 5 
        f = [2, 3, 7, 7];
        A = [1, 1, -2, -5;-1, 2, 1, 4];
        b = [2, 3]';
        e = [1, 1];
        lb = zeros(4,1); 
        ub = [30 100 20 1]';  
        xint = 'CICI';
        prob = optiprob('f',f,'mix',A,b,e,'bounds',lb,ub,'int',xint);            
        sol = [0;2;0;0];
        fmin = 6;
        
    case 6 
        f = [1, 2, 3, 7, 8, 8];
        A = [5, -3, 2, -3, -1, 2; -1, 0, 2, 1, 3, -3;1, 2, -1, 0, 5, -1];
        b = [-5, -1, 3]';
        e = [1, 1, 1];
        lb = zeros(6,1);
        ub = 10*ones(6,1); 
        xint = 'IIIIII';
        prob = optiprob('f',f,'mix',A,b,e,'bounds',lb,ub,'int',xint);            
        sol = [1;1;0;0;0;0];
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
        xint = [1 3];
        prob = optiprob('f',f,'mix',A,b,e,'bounds',lb,ub,'int',xint);            
        sol = [4;0;0;0.2;zeros(40,1)];
        fmin = 4.0000000000000044;
        
    case 8 
        f = -[8, 15];
        A = [10, 21;2, 1];
        b = [156, 22]';
        e = [-1, -1]; 
        xint = 'IC';
        prob = optiprob('f',f,'mix',A,b,e,'int',xint);            
        sol = [9;3.14285714285714];
        fmin = -119.14285714285714;
        
    case 9 
        f = -[3, 13];
        A = [2, 9;11, -8];
        b = [40, 82]';
        e = [-1, -1]; 
        xint = 'IC';
        prob = optiprob('f',f,'mix',A,b,e,'int',xint);            
        sol = [9;2.44444444444444];
        fmin = -58.77777777777778;
        
    case 10 
        f = -[592, 381, 273, 55, 48, 37, 23];
        A = [3534, 2356, 1767, 589, 528, 451, 304];
        b = 119567;
        e = -1;
        lb = zeros(7,1); 
        ub = [100 50 33 20 77 44 20]';
        xint = 'IIIIIII';
        prob = optiprob('f',f,'mix',A,b,e,'bounds',lb,ub,'int',xint);            
        sol = [32;2;1;0;0;0;0];
        fmin = -19979;        
        
    otherwise
        error('Problem not available or not implemented yet');
end