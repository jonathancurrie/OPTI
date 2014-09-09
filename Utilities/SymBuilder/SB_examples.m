%% Few optimization test problems for SymBuilder
clear all
clc
%Unfinished code, but sort of working! Requires the Symbolic Toolbox.

%% LP1 [fval = -31.4]
clc
B = SymBuilder();

B.AddObj('-6*x1 -5*x2');
B.AddCon('x1 + 4*x2 <= 16');
B.AddCon('6*x1 + 4*x2 <= 28');
B.AddCon('2*x1 - 5*x2 <= 6');
B.AddBound('0 <= x <= 10');

Build(B)

[x,fval,ef,info] = Solve(B)

%% MILP1 [fval = -15]
clc
B = SymBuilder();

B.AddObj('3*x1 -7*x2 -12*x3');
B.AddCon('-3*x1 + 6*x2 + 8*x3 <= 12');
B.AddCon('6*x1 - 3*x2 + 7*x3 <= 8');
B.AddCon('-6*x1 + 3*x2 + 3*x3 <= 5');
B.AddInteger('x = I');

Build(B)

[x,fval,ef,info] = Solve(B)

%% QP1 [fval = -6.4138]
clc
B = SymBuilder();

B.AddObj('0.5*x1^2 -1*x1*x2 + x2^2 - 2*x1 - 6*x2');
B.AddCon('x1 + x2 <= 2');
B.AddCon('-x1 + 2*x2 <= 2');
B.AddCon('2*x1 + x2 <= 3');
B.AddCon('x1 + 1.5*x2 = 2');
B.AddBound('0 <= x <= 10');

Build(B)

[x,fval,ef,info] = Solve(B)

%% MIQP1 [fval = -2.75]
clc
B = SymBuilder();

B.AddObj('0.5*x1^2 + 0.5*x2^2 + 0.5*x3^2 - 2*x1 - 3*x2 - x3');
B.AddCon('x1 + x2 + x3 <= 1');
B.AddCon('3*x1 - 2*x2 - 3*x3 <= 1');
B.AddCon('x1 - 3*x2 + 2*x2 <= 1');
B.AddInteger('x2 = I');

Build(B)

[x,fval,ef,info] = Solve(B)

%% UNO [fval = 0]
clc
B = SymBuilder();

B.AddObj('(1-x1)^2 + 100 * (x2-x1^2)^2');

Build(B)

[x,fval,ef,info] = Solve(B,[0;0],symbset('cbmode','cppad','display','iter'))

%% NLP1 [fval = 0]
clc
B = SymBuilder();

B.AddObj('(x1-x2)^2 + (x2+x3-2)^2 + (x4-1)^2 + (x5-1)^2');
B.AddCon('x1 + 3*x2 = 4');
B.AddCon('x3 + x4 - 2*x5 = 0');
B.AddCon('x2 - x5 = 0');

Draft(B)

x0 = [ 2.5 0.5 2 -1 0.5 ];
[x,fval,ef,info] = Solve(B,x0,symbset('cbmode','cppad','display','iter'))


%% MINLP1 [fval = -0.72007]
clc
clear
B = SymBuilder(false);

B.AddObj('sin(pi*x1/12)*cos(pi*x2/16)');
B.AddCon('-x1 + 2.5*x2 <= 1');
B.AddCon('x1 + 2.5*x2 <= -15');
B.AddInteger('x = I');

Build(B)
Opt = GetOPTI(B,symbset('cbmode','mcode'));
Opt = opti(Opt,optiset('derivCheck','on'))
x0 = [0;0];
[x,fval,ef,info] = solve(Opt,x0)

%% Changing Solver 
clc
B = SymBuilder();

B.AddObj('-6*x1 -5*x2');
B.AddCon('x1 + 4*x2 <= 16');
B.AddCon('6*x1 + 4*x2 <= 28');
B.AddCon('2*x1 - 5*x2 <= 6');
B.AddBound('0 <= x <= 10');

Build(B)
opts = symbset('solver','scip');

[x,fval,ef,info] = Solve(B,[],opts)

%% Approximate 2nd Derivatives
clc
B = SymBuilder();

B.AddObj('(x1-x2)^2 + (x2+x3-2)^2 + (x4-1)^2 + (x5-1)^2');
B.AddCon('x1 + 3*x2 = 4');
B.AddCon('x3 + x4 - 2*x5 = 0');
B.AddCon('x2 - x5 = 0');
B.AddCon('x2^2 + x5 = 4');
B.AddInteger('x3 = I');

Build(B)

x0 = [ 2.5 0.5 2 -1 0.5 ];
opts = symbset('use2ndDerivs','no','cbmode','mcode');
[x,fval,ef,info] = Solve(B,x0,opts)


%% NLP1 C CODE
clc
B = SymBuilder();

B.AddObj('(x1-x2)^3 + (x2+x3-2)^2 + (x4-1)^2 + (x5-1)^2');
B.AddCon('x1 + 3*x2 = 4');
B.AddCon('x3 + x4 - 2*x5 = 0');
B.AddCon('x2 - x5^2.2 = 0');

Build(B)


opts = symbset('use2ndDerivs','yes','cbmode','auto');
x0 = [ 2.5 0.5 2 -1 0.5 ];
sig = 1.1; lambda = [1;2;3];
prob = GetNLProb(B,opts);
% Opt = GetOPTI(B,opts);
% oopts = optiset('derivCheck','on');
% Opt = opti(Opt,oopts)

prob.f(x0)
full(prob.nljac(x0))
full(prob.H(x0,sig,lambda))


%% Mixed Var Name NLP C CODE
clc
B = SymBuilder();

B.AddObj('(btb-m1)^3 + (m1+x3-2)^2 + (x4-1)^2 + (x5-1)^2');
B.AddCon('btb + 3*m1 = 4');
B.AddCon('x3 + x4 - 2*x5 = 0');
B.AddCon('m1 - x5^2.2 = 0');

Build(B)


opts = symbset('use2ndDerivs','yes','cbmode','mcode');
x0 = [ 2.5 0.5 2 -1 0.5 ];
sig = 1.1; lambda = [1;2;3];
prob = GetNLProb(B,opts);
% Opt = GetOPTI(B,opts);
% oopts = optiset('derivCheck','on');
% Opt = opti(Opt,oopts)

prob.f(x0)
% symb_cb('grad',x0)
full(prob.nljac(x0))
% full(symb_cb('jac',x0))
full(prob.H(x0,sig,lambda))
% full(symb_cb('hess',x0,sig,lambda))


%% ALTERNATIVE SETUP EXAMPLES
%% NLP Example
clc
%Create a Symbolic Builder Object
B = SymBuilder();

%Add Optimization Problem
B + 'min (x1-x2)^2 + (x2+x3-2)^2 + (x4-1)^2 + (x5-1)^2';
B + 'x1 + 3*x2 = 4';
B + 'x3 + x4 - 2*x5 = 0';
B + 'x2 - x5 = 0';
%Construct The Problem
Draft(B)
%Solve it using IPOPT
x0 = [ 2.5 0.5 2 -1 0.5 ];
[x,fval,ef,info] = Solve(B,x0,symbset('use2ndDerivs','no'))


%% LP1 [fval = -31.4]
clc
B = SymBuilder();

B + 'min -6*x1 -5*x2';
B + 'x1 + 4*x2 <= 16';
B + '6*x1 + 4*x2 <= 28';
B + '2*x1 - 5*x2 <= 6';
B + '0 <= x <= 10';

Build(B)

[x,fval,ef,info] = Solve(B)

%% LP1 [fval = -31.4]
clc
B = SymBuilder();
A = [1 4; 6 4; 2 -5];
rl = -Inf(3,1); ru = [16;28;6];

B + 'min -6*x1 -5*x2';
B.AddLinCon(sym('[x1;x2]'),A,rl,ru);
B + '0 <= x <= 10';

Build(B)

[x,fval,ef,info] = Solve(B)

%% MILP1 [fval = -15]
clc
B = SymBuilder();

B + 'min 3*x1 -7*x2 -12*x3';
B + '-3*x1 + 6*x2 + 8*x3 <= 12';
B + '6*x1 - 3*x2 + 7*x3 <= 8';
B + '-6*x1 + 3*x2 + 3*x3 <= 5';
B + 'x = I';

Build(B)

[x,fval,ef,info] = Solve(B)

%% QP1 [fval = -6.4138]
clc
B = SymBuilder();

B + 'min 0.5*x1^2 -1*x1*x2 + x2^2 - 2*x1 - 6*x2';
B + 'x1 + x2 <= 2';
B + '-x1 + 2*x2 <= 2';
B + '2*x1 + x2 <= 3';
B + 'x1 + 1.5*x2 = 2';
B + '0 <= x <= 10';
Build(B)

[x,fval,ef,info] = Solve(B)


%% QP1 with Symbolic Variables
clc
B = SymBuilder(false);

obj = @(x) (x(1)-x(2))^2 + (x(2)+x(3)-2)^2 + (x(4)-1)^2 + (x(5)-1)^2;
nlcon = @(x) [x(1) + 3*x(2); x(3) + x(4) - 2*x(5); x(2) - x(5)];
cl = [4;0;0];
cu = cl;

%Evaluate with Symbolic Variables
x = sym('x',[5,1]);
sym_obj = obj(x);
sym_nlcon = nlcon(x);

%Add to Symbuilder
B.AddObj(sym_obj);
B.AddCon(sym_nlcon,cl,cu);

Build(B);

%Solve
x0 = [ 2.5 0.5 2 -1 0.5 ];
[x,fval,ef,info] = Solve(B,x0,symbset('cbmode','cppad'))

