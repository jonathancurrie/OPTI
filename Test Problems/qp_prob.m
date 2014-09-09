function [prob,sol,fmin] = qp_prob(varargin)
%QP_PROB  Return an OPTI QP 
%
%   prob = qp_prob(no) return a pre-built optiprob of a saved QP.
%
%   [prob,sol,fmin] = qp_prob(no) returns the optimum solution and 
%   function eval at the optimum
%
%   no = qp_prob() returns the number of problems available for testing.

%   (C) 2011 Jonathan Currie (I2C2)

%Check if just returning no problems
if(nargin < 1)
    prob = 3; sol = []; fmin = [];
    return;
else
    no = varargin{1};
end          

%Big switch yard
switch(no)
    case 1 
        H = eye(3);
        f = -[2 3 1]';
        A = [1 1 1;3 -2 -3; 1 -3 2]; 
        b = [1;1;1];
        prob = optiprob('H',H,'f',f,'ineq',A,b);            
        sol = [1/3;4/3;-2/3];
        fmin = -2.83333333301227;
        
    case 2 
        H = [1 -1; -1 2];
        f = -[2 6]';
        A = [1 1; -1 2; 2 1];
        b = [2; 2; 3]; 
        lb = [0;0];
        prob = optiprob('H',H,'f',f,'ineq',A,b,'lb',lb);             
        sol = [1/3;4/3];
        fmin = -8.22222220552525;
        
    case 3 
        H = [1 -1; -1 2];
        f = -[2 6]';
        A = [1 1; -1 2; 2 1];
        b = [2; 2; 3];
        Aeq = [1 1.5];
        beq = 2;
        lb = [0;0];
        ub = [10;10];
        prob = optiprob('H',H,'f',f,'ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub);           
        sol = [0.34482763667623;1.10344824221585];
        fmin = -6.41379310344827;               
        
    otherwise
        error('Problem not available or not implemented yet');
end