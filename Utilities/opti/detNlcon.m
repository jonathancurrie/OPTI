function [nlcon,nlrhs,nle] = detNlcon(nonlcon,x0)
%DETNLCON  Determine OPTI format Nonlinear Constraints

%   Copyright (C) 2011 Jonathan Currie (I2C2)

%Sort out nonlinear constraints (unfortunately rather inefficient)
if(isempty(nonlcon))
    nlcon = []; nlrhs = []; nle = [];
    return;
end

[in,eq] = nonlcon(x0);
nineq = length(in);
neq = length(eq);

nlcon = @(x) conval(nonlcon,x);
nlrhs = zeros(nineq+neq,1);
nle = [-1*ones(nineq,1);zeros(neq,1)];


%Dummy function
function c = conval(nonlcon,x)
[in,eq] = nonlcon(x);
c = [in;eq];

