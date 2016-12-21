function [f,ind] = func2sym(fun,vars)
%FUNC2SYM  Convert a function handle to symbolic variable equation
%
%   f = func2sym(fun) converts the function handle into a string, then
%   replaces any x(1) to x1, x(2) to x2, etc, and returns it as a symbolic
%   equation.
%
%   [f,ind] = func2sym(fun) returns the indices of the symbolic variable x
%   used in the equation.

%   Copyright (C) 2013 Jonathan Currie (IPL)

if(nargin < 2), vars = {'x'}; end
if(~iscell(vars)), vars = {vars}; end

%Create Outputs
ind = cell(size(vars));

%Convert function to string
f = func2str(fun);

%Convert each variable into 'symbolic var' (e.g. x(1) to x1, z(1) to z1)
for i = 1:length(vars)
    [f,ind{i}] = buildSym(f,vars{i},i~=1);
end

%Convert to symbolic string
f = sym(f);

%Substitute any constants from original workspace into symbolic expression
fdata = functions(fun);
if(isfield(fdata,'workspace') && ~isempty(fdata.workspace))
    wrkspc = fdata.workspace{1};
    f = subs(f,fieldnames(wrkspc),struct2cell(wrkspc));
end


function [f,ind] = buildSym(str,var,skipStart)
if(nargin < 3), skipStart = false; end

k = strfind(str,[var '(']);
%Check for exp(, remove as required
k2 = strfind(str,'exp('); k2 = k2 + 2;
k = setdiff(k,k2);
len = length(k);
ind = zeros(len,1);

%Determine X Indicies in order
for i = 1:len
    ind(i) = readVar(str,k(i));
end

%Start building symbolic string
if(skipStart)
    start = 0;
else
    start = findBrac(str,1); %get first bracket @(x)
end
f = '';
for i = 1:len
    f = [f str(start+1:k(i)) num2str(ind(i))]; %#ok<AGROW>
    start = findBrac(str,k(i));
end
%Add remainder of string
f = [f str(start+1:end)];


function ind = readVar(str,k)
start = k+2; %skip x(
stop = findBrac(str,start)-1;
ind = str2double(str(start:stop));

function ind = findBrac(str,k)
for i = 0:length(str)
    if(str(k+i) == ')')
        ind = k+i;
        break
    end
end