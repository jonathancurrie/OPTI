function [prob,sol,fmin] = bilp_prob(varargin)
%BILP_PROB  Return an OPTI BILP 
%
%   prob = bilp_prob(no) return a pre-built optiprob of a saved BILP.
%
%   [prob,sol,fmin] = bilp_prob(no) returns the optimum solution and 
%   function eval at the optimum
%
%   no = bilp_prob() returns the number of problems available for testing.

%   (C) 2011 Jonathan Currie (IPL)

%Check if just returning no problems
if(nargin < 1)
    prob = 2; sol = []; fmin = [];
    return;
else
    no = varargin{1};
end          

%Big switch yard
switch(no)
    case 1 
        f = -[6 5]';
        A = [-3,5; 6,4; 3, -5; -6, -4]; 
        b = [6;9;1;3];
        xint = 'BB';
        prob = optiprob('f',f,'ineq',A,b,'int',xint);            
        sol = [0;1];
        fmin = -5;
        
    case 2 
        f = -[9 5 6 4]';
        A = [6 3 5 2; 0 0 1 1; -1 0 1 0; 0 -1 0 1];
        b = [9; 1; 0; 0];
        xint = 'BBBB';
        prob = optiprob('f',f,'ineq',A,b,'int',xint);              
        sol = [1;1;0;0];
        fmin = -14;                     
        
    otherwise
        error('Problem not available or not implemented yet');
end