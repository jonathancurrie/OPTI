function [prob,sol,fmin] = nlp_prob(varargin)
%NLP_PROB  Return an OPTI NLP 
%
%   prob = nlp_prob(no) return a pre-built optiprob of a saved NLP.
%
%   [prob,sol,fmin] = nlp_prob(no) returns the optimum solution and 
%   function eval at the optimum
%
%   no = nlp_prob() returns the number of problems available for testing.
%
%   NOTE as implemented problems 1-50 return a Hock-Schittkowski
%   problem (see nlp_HS), while other numbers are other problems.

%   (C) 2011 Jonathan Currie (I2C2)

%Check if just returning no problems
if(nargin < 1)
    prob = 50; sol = []; fmin = [];
    return;
else
    no = varargin{1};
end       

%Big switch yard
switch(no)
    case num2cell(1:50) %HS Problems
        [prob,sol,fmin] = nlp_HS(no);
                            
        
    otherwise
        error('Problem not available or not implemented yet');
end