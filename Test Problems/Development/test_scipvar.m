%% Testing SCIPVAR class

clear all
clc

%% T1
clear all
clc

x = scipvar(1);
t1 = 2*x

%% T2
clear all
clc

x = scipvar(1);
t1 = x*2

%% T3
clear all
clc

x = scipvar(1);
t1 = 2+x*3

%% T4
clear all
clc

x = scipvar(1);
t1 = x*3+2

%% T5
clear all
clc

x = scipvar(1);
t1 = 2-x*3

%% T5
clear all
clc

x = scipvar(1);
t1 = x*3 - 2

%% T6
clear all
clc

x = scipvar(1);
t1 = 2/x*3

%% T7
clear all
clc

x = scipvar(1);
t1 = 2/(x*3)

%% T8
clear all
clc

x = scipvar(1);
t1 = (2/x)*3

%% T9
clear all
clc

x = scipvar(1);
t1 = x*3 / (2*x)

%% T10
clear all
clc

x = scipvar(1);
t1 = (1-x) / (1+x)*2

%% T11
clear all
clc

x = scipvar(1);
t1 = (1-x) / ((1+x)*2)

%% T12
clear all
clc

x = scipvar(1);
t1 = 1-(x+2)

%% T13
clear all
clc

x = scipvar(2,1)

t1 = 3*x(1)*x(2) + 2*x(1) - 3*x(2)/4


%%
clear all
clc

x = scipvar(2,1)

t2 = [3*x(1)*x(2) + 2*x(1) - 3*x(2)/4
      2*x(1)*x(2)]

%% Vectorized Testing

%% Vection addition1 
%double + scipvar
clc
clear all
a = 2;
b = [1 3 5; 2 4 6];
x = scipvar(2,3);

%2+x1
y2 = a+x(1)

%2+[x..]
y = a+x

%[1..]+x1
y3 = b+x(1)

%[1..]+[x..]
z = b+x

z2 = 2 - b.*x

%% Vector addition2
%scipvar + double
clc
clear all
a = 2;
b = [1 3 5; 2 4 6];
x = scipvar(2,3);

%x1+2
x(1)+a

%[x..]+2
y = x+a

%x1+[1..]
y2 = x(1)+b

%[x..]+[1..]
z = x+b

%% Vector addition3
%scipvar + scipvar
clc
clear all
x = scipvar(2,3);

%x1+x1
y = x(1)+x(1)

%x1+[x..]
y2 = x(2)+x

%[x..]+x1
z = x+x(1)

%[x..]+[x..]
z2 = x+x

z3 = 2*x - x

%% Vectorized Functions
%scipvar + scipvar
clc
clear all
x = scipvar(2,3);

y = log(x)

%% Dot product
clc
clear all
n = 10;
x = scipvar(n,1);

y = dot(x,x);
x0 = rand(n,1);

opti_scipnltest(@(x) dot(x,x),x0)

%% Matrix vector0
%double x scipvar
clc
clear all
X = scipvar(2,2);
x = scipvar(2,1);
h = [1;2];
H = [1 2; 3 4];
H2 = [1 2; 3 4; 5 6];
X2 = scipvar(2,3);

%H*x(1)
z = H*x(1)

%H*x
z2 = H*x

%H*X
z3 = H*X

%Uneven sizes
z4 = H2*X2

%Inner product
z5 = h'*x

%Outer product
z6 = h*x'

%% Matrix vector1
%double x scipvar
clc
clear all
X = scipvar(2,2);
x = scipvar(2,1);
h = [1;2];
H = [1 2; 3 4];
H2 = [1 2; 3 4; 5 6];
X2 = scipvar(2,3);

%H*x(1)
z = H*x(1)

%H*x
z2 = H*x

%H*X
z3 = H*X

%Uneven sizes
z4 = H2*X2

%Inner product
z5 = h'*x

%Outer product
z6 = h*x'

%% Matrix vector2
%scipvar x double
clc
clear
X = scipvar(2,2);
x = scipvar(2,1);
h = [1;2];
H = [1 2; 3 4];
H2 = [1 2; 3 4; 5 6; 7 8];
X2 = scipvar(2,4);

% x(1)*H
x(1)*H

%x'*H
x'*H

%X*H
X*H

%Uneven sizes
X2*H2

%Inner product
x'*h

%Outer product
x*h'

%% Matrix vector3
%scipvar x scipvar
clc
clear
X = scipvar(2,2);
x = scipvar(2,1);
X2 = scipvar(2,4);

% x(1)*x(1)
x(1)*x(1)

%x'*X
x'*X

%X*H
X*X

%Uneven sizes
X2*X2'

%Inner product
x'*x

%Outer product
x*x'


%% Loop Checking
clc
% clear all
b = 1:3;
x = scipvar(3,1);
j = x(1)*x(2);
for i = 1:3
    j = j + b(i)*x(i);
end
j

%% Test1 2D
clear all
clc
a = 100;
x = scipvar(2,1);

fun = @(x) a*(x(2) - x(1)^2)^2 + (1 - x(1))^2
fun = @(x) 1/(27*sqrt(3)) * ((x(1) - 3)^2 - 9)*x(2)^3;
fun = @(x) 0.01*(x(1)-1)^2 + (x(2)-x(1)^2)^2;

b = fun(x)

x0 = [-1;2];
fun(x0)
% eval(b,x0)
opti_scipnltest(fun,x0)

%% Test2 3D
clear all
clc
a = 100;
x = scipvar(3,1);

fun = @(x) (x(1)-1)*(x(1)-2)*(x(1)-3)+x(3);
fun = @(x) 9 - 8*x(1) - 6*x(2) - 4*x(3) + 2*x(1)^2 + 2*x(2)^2 + x(3)^2 + 2*x(1)*x(2) + 2*x(1)*x(3);

b = fun(x)

x0 = [-1;2;4];
fun(x0)
% eval(b,x0)
opti_scipnltest(fun,x0)

%% Test3 4D
clear all
clc
a = 100;
x = scipvar(4,1);

fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2 + 90*(x(4)-x(3)^2)^2 + (1-x(3))^2 + 10.1*((x(2)-1)^2 + (x(4)-1)^2) + 19.8*(x(2)-1)*(x(4)-1);
fun = @(x) -x(1)*x(2)*x(3)*x(4);
fun = @(x) 2 - x(1)*x(2)*x(3);
fun = @(x) x(1)^2 + x(2)^2 + 2*x(3)^2 + x(4)^2 - 5*x(1) - 5*x(2) - 21*x(3) + 7*x(4);

b = fun(x)

x0 = [-1;2;4;-5];
fun(x0)
% eval(b,x0)
opti_scipnltest(fun,x0)

%% Test4 5D
clear all
clc
a = 100;
x = scipvar(5,1);

fun = @(x) (x(1)-x(2))^2 + (x(2)-x(3))^3 + (x(3)-x(4))^4 + (x(4)-x(5))^4;

b = fun(x)

x0 = [-1;2;4;-5;7];
fun(x0)
% eval(b,x0)
opti_scipnltest(fun,x0)

%% PROD test
clc
clear all
x = scipvar(3,1);
x1 = scipvar(1,3);
x2 = scipvar(3);
x0 = [1:3;4:6;7:9];

prod(x)
prod(x1)
y = prod(x2,1)
y2 = prod(x2,2)
prod(y)

x0
prod(x0,1)
prod(x0,2)

opti_scipnltest(@(x) prod(prod(x)),x0)

%% SUM test
clc
clear all
x = scipvar(3,1);
x1 = scipvar(1,3);
x2 = scipvar(3);
x0 = [1:3;4:6;7:9];

z = sum(x)
z2 = sum(x1)
y = sum(x2,1)
y2 = sum(x2,2)
sum(y)

x0
sum(x0,1)
sum(x0,2)

opti_scipnltest(@(x) sum(sum(x)),x0)

%% NORM test
%note 2norm for vectors, Frobenius norm for matrices
clc
clear all
x = scipvar(3,1);
X = scipvar(3);

norm(x)
norm(X)

%% MPOWER integer
clc
clear all
x = scipvar(5);
x1 = scipvar(2,1);

x^2

x1.^2

%% Function checking
clc
clear all
a = 2;
b = (1:3)';
x = scipvar(3,1);
x0 = ones(3,1);

fun = @(x) prod(b.*x);

opti_scipnltest(fun,x0)

%% Special Power 0
clear
clc
x = scipvar(2,1);

f = @(x) exp((3-x(2)).*log(x(1)));
% f = @(x) (3-x(2)).*log(x(1));
% f = @(x) x(2) - exp(x(1))
y = f(x)

x0 = [0.5;1.5];
f(x0)

opti_scipnltest(f,x0)

%% Special Power
clear
clc
x = scipvar(2,1);

f = @(x) x(1)^(3-x(2));
y = f(x)

x0 = [0.5;14]; 
f(x0)
% eval(y,x0)
opti_scipnltest(f,x0)

%% Special Power Vec
clear
clc
x = scipvar(2,1);

f = @(x) sum(2*x.^(4-x(2)))
y = f(x)

x0 = [0.5;14]; 
f(x0)
% eval(y,x0)
opti_scipnltest(f,x0)

%% Diag
clc
clear all
x = scipvar(2,1);
X = scipvar(2,2);
X2 = scipvar(2,4);
X3 = scipvar(3,3);

H = [1 3; 2 4]
H2 = [1 4 7; 2 5 8; 3 6 9];
diag(H)
diag(H,1)
diag(H,-1)
diag(H,2)
diag(H2,1)
diag(H2,2)

diag(X)
diag(X,1)
diag(X,-1)
diag(X2)
diag(X3)
diag(X3,1)
diag(X3,2)
diag(X3,-2)

diag(x)

%% TRIU Test
clc
clear
x = scipvar(5);

y = triu(x)

y2 = triu(x,1);

y3 = triu(x,-1);

%% TRIL Test
clc
clear
x = scipvar(5);

y = tril(x)

y2 = tril(x,1);

y3 = tril(x,-1);

%% Colon Test
clc
clear
x = scipvar(3);
x2 = scipvar(1,3)

x(:)

x2(:)

%% Constraint Test
clc
clear
a = 100;
x = scipvar(5,1);

nlcon = @(x) [x(1) + 2*x(2) + 3*x(3) - 6;
              x(2) + 2*x(3) + 3*x(4) - 6;
              x(3) + 2*x(4) + 3*x(5) - 6];

b = nlcon(x)





%% Constraint Test 2
clc
clear 
a = 100;
x = scipvar(5,1);

c(3) = x(1) + 2*x(2) + 3*x(3) - 6;
c(2) = x(2) + 2*x(3) + 3*x(4) - 6;
c(1) = x(3) + 2*x(4) + 3*x(5) - 6;
c

y = c*x(1:3)

%% Constraint Test 3
clc
clear all
a = 100;
x = scipvar(5,1);

c(1,1) = x(1) + 2*x(2) + 3*x(3) - 6


%% Sparse Testing
clc
clear all
x = scipvar(3,1);
A = speye(3);

y = x'*A*x

%% Nick S. Transport Problem
clc
clear

%Problem
a = [350;600]; %capacity of plants
b = [325;300;275]; %demand at market
d = [2.5 1.7 1.8; %distance table
     2.5 1.8 1.4];

fc = 90; %cost per thousand miles
c = fc*d/1000; %transport cost per case

%Objective
fun = @(x) sum(sum(c.*x));

%Constraints
nlcon = @(x) [sum(x,2);   
              sum(x,1)'];  
cl = [-Inf(2,1); b];
cu = [a; Inf(3,1)];  
%Bounds
lb = zeros(6,1);
%Starting guess (sets X size)
x0 = zeros(2,3);

%Solve
opts = optiset('display','iter');
[x2,f2,e2,i2] = opti_scipnl(fun,[],[],[],lb,[],nlcon,cl,cu,[],x0,opts) 
