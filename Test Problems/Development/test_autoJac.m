%% Testing Will's adiff
clc
clear
%Array of functions to check
funs = {'abs','acos','asin','atan','cos','cosh','exp','log','log10','sin','sinh','sqrt','tan','tanh'};

%Check against both numerical and symbolic gradients
n = length(funs);
ajac = zeros(n,1);
njac = zeros(n,1);
sjac = zeros(n,1);

x0 = abs(0.5+0.2*randn);

%Test each one
for i = 1:n
    fun = str2func(['@(x)' funs{i} '(x)'])
    ajac(i) = autoJac(fun,x0);
    njac(i) = mklJac(fun,x0);
    s = symJac(fun);
    sjac(i) = s(x0);
end

%Display Results
raw = [ajac njac sjac]
acc = [abs(sjac-ajac) abs(sjac-njac)]

%% Index Checking

%% Vection addition1 
%double + adiff
clc
clear 
a = 2;
b = [1 3 5; 2 4 6];
x = adiff(2);

%2+x1
y2 = a+x(1)

%2+[x..]
y = a+x

%[1..]+x1
y = b+x(1)

%[1..]+[x..]
z = b+x

%% Vector addition2
%adiff + double
clc
clear
a = 2;
b = [1 3 5; 2 4 6];
x = adiff(2);

%x1+2
x(1)+a

%[x..]+2
y = x+a

%x1+[1..]
x(1)+b

%[x..]+[1..]
z = x+b

%% Vector addition3
%adiff + adiff
clc
clear all
x = adiff(2);

%x1+x1
y = x(1)+x(1)

%x1+[x..]
y2 = x(1)+x

%[x..]+x1
z = x+x(1)

%[x..]+[x..]
z2 = x+x

%% Loop Checking
clc
% clear all
b = 1:3;
x = adiff(1:3);
j = x(1)*x(2);
for i = 1:3
    j = j + b(i)*x(i);
end
j

%% Test1 2D
clear all
clc
x0 = [-1;2];
x = adiff(x0);

fun = @(x) 1/(27*sqrt(3)) * ((x(1) - 3)^2 - 9)*x(2)^3;

b = fun(x)

fun(x0)
[f,j] = adiffget(b)

s = symJac(fun);
j2 = s(x0)

%% Test2 3D
clear all
clc
x0 = [-1;2;4];
x = adiff(x0);

fun = @(x) 9 - 8*x(1) - 6*x(2) - 4*x(3) + 2*x(1)^2 + 2*x(2)^2 + x(3)^2 + 2*x(1)*x(2) + 2*x(1)*x(3);

b = fun(x)


fun(x0)
[f,j] = adiffget(b)

s = symJac(fun);
j2 = s(x0)

%% Test3 4D
clear all
clc
x0 = [-1;2;4;-5];
x = adiff(x0);

fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2 + 90*(x(4)-x(3)^2)^2 + (1-x(3))^2 + 10.1*((x(2)-1)^2 + (x(4)-1)^2) + 19.8*(x(2)-1)*(x(4)-1);

b = fun(x)

fun(x0)
[f,j] = adiffget(b)

s = symJac(fun);
j2 = s(x0)

%% Test4 5D
clear all
clc
a = 100;
x0 = [-1;2;4;-5;7];
x = adiff(x0);

fun = @(x) (x(1)-x(2))^2 + (x(2)-x(3))^3 + (x(3)-x(4))^4 + (x(4)-x(5))^4;

fun(x0)
b = fun(x)

[f,j] = adiffget(b)

s = symJac(fun);
j2 = s(x0)

%% Subsasgn, Empty LHS
clear all
clc
x = adiff(1:4);
X = adiff(2);
a = [1;2;3;4];

%Simple indexing
y = x(1)*x(2)

%Simple assign
y2(2) = x(3)*2*x(4)
%Simple 2D assign
y22(2) = x(3)*2*x(4)

% %Multi index, indexing
z = x(1:3)

%Multi index, assign
y3(1:2) = [x(1);x(3)]

%Multi index, assign, longer vec
y4(2:3) = [x(3);x(4)]

%Logical Index
y5(a>2) = x(a>2)

%% Subsasgn, Existing LHS
clear
clc
x = adiff(1:4);
a = [1;2;3;4];

%Add element to vector
y = x(1)
y(2) = x(2)
y(3) = x(3)
y(5) = x(4)

%Add element to indexed vector
y2(1) = x(1)

%Add element to existing vector
% y2(:,2) = x(1:2)
% y2(:,3) = x;


%% Subsagn, Element Numbers
clear
clc
x = adiff(2);

% X([8 9]) = [x(1); x(2)]

% X(:,1:2) = barvec(2)


%% Multiline Equation + multi index
clear
clc
x = adiff(1:4);

f = 3*x.^2 - x(1)
g = sum(f)
z1(1:4) = f
% z1(5) = g
% y = z1(5)*2
% 
% %Multi Assign
% c(1,1) = 2*x(1);
% c(2:3,1) = 3*x(2:3).^2;
% c(4,1) = 0.1*x(4)


%% Multi index subs
clear
clc
x = adiff(1:3);

y = sum(x(1:2).*x(2:3))

% y.vecstr

%% Multi2
clear all
clc

b = [1:4]'
x = adiff(1:4);

f = @(x) b(1:2).*x(1:2) + b(3:end).*x(3:end).^2
y = f(x)


x0 = ones(4,1); 
f(x0)

%% Multi assign
clear
clc

x = adiff(1:4);

f([2 3]) = 2*x([3 2])
f([1 4]) = 5*x([1 4])

%% Constraint Test
clc
clear
x = adiff(1:5);

nlcon = @(x) [x(1) + 2*x(2) + 3*x(3) - 6;
              x(2) + 2*x(3) + 3*x(4) - 6;
              x(3) + 2*x(4) + 3*x(5) - 6];

b = nlcon(x)

[f,j] = adiffget(b)

s = symJac(nlcon);
j2 = s(ones(5,1))

%% Constraint Test 2
clc
clear 
a = 100;
x = adiff(1:5);

c(3) = x(1) + 2*x(2) + 3*x(3) - 6;
c(2) = x(2) + 2*x(3) + 3*x(4) - 6;
c(1) = x(3) + 2*x(4) + 3*x(5) - 6;
c

% y = c*x(1:3)

[f,j] = adiffget(c)

%% DIW

clc
clear all

x = adiff([1 1])

length(x)
size(x)
ndims(x)

n = length(x); 
j = @(x) sum((1-x).^2) + sum(100*(x(2:n) - x(1:n-1).^2).^2)


[f,jac] = adiffget(j(x))

mklJac(j,ones(2,1))
