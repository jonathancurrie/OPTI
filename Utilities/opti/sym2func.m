function f = sym2func(fun,args)
%SYM2FUNC  Convert a symbolic variable equation to function handle
%
%   f = sym2func(fun) converts the symbolic equation into a string, then
%   replaces any x1 to x(1), x2 to x(2), etc, and returns it as a function
%   handle suitable for standard evaluation.

%   Copyright (C) 2013 Jonathan Currie (IPL)

if(nargin < 2), args = 'x'; end

%Extract symbolic variables in the expression
sv = symvar(fun);
%For each symbolic var, convert x1 to x(1) as cell of strings
cv = cell(size(sv));
for i = 1:length(sv)
    cv{i} = symVar2MVar(sv(i));
end
%Substitute our indexed variables into the symbolic expression
sfun = subs(fun,sv,cv);

%Convert to string
str = char(sfun);
%Remove 'matrix' if it is found
k = strfind(str,'matrix');
if(k)
    str = str(k+7:end-1);
end
%Ensure we have right size
[r,c] = size(fun);
if(c==1 && r > 1)
    str = regexprep(str,',',';'); %force column
elseif(c>1 && r>1)
    str = regexprep(str,'],','];'); %keep matrix
end

if(iscell(args))
    a = args{1};
    for i = 2:length(args)
        a = sprintf('%s,%s',a,args{i});
    end
    args = a;
end

%Build function handle
f = eval(['@(' args ') ' str]);


function str = symVar2MVar(svar)
s = char(svar);
for i = 2:length(s)
    if(~isalpha(s(i)))
        str = [s(1:i-1) '(' s(i:end) ')'];
        return;
    end
end
%No index?, oh well... could be for e.g. just t
str = s;


function is = isalpha(str)
if((str>='a' && str<='z') || (str>='A' && str<='Z'))
    is = true;
else
    is = false;
end
