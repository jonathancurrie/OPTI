%%
clear
load lp_test_results

%%
n = 1
lp = lp_tprob{n}
lp_sval(n)

x = cplexlp(lp.f,lp.A,lp.b,lp.Aeq,lp.beq,lp.lb,lp.ub,lp.int)

[A,rl,ru] = gen2row(lp.A,lp.b,lp.Aeq,lp.beq);

%%
clc
[x,f,e,i] = opti_cplex(lp.H,lp.f,A,rl,ru,lp.lb,lp.ub,lp.int)



%%
clc
%Run the following line to see what solvers are available to OPTI:
checkSolver();

%Alternatively you can find the best available solver for a given problem 
%type:
lp = checkSolver('best_lp')

%Or see all available solvers for a given problem type
checkSolver('lp')
% for i = 1:1e3
    f = -[6 5]';
    A = ([1,4; 6,4; 2, -5]); 
    b = [16;28;6];    
    lb = [0;0];
    ub = [10;10];

    cp = Cplex('opti');
    cp.Model.obj = f;
    cp.Model.A = A;
    cp.Model.lhs = -Inf(size(b));
    cp.Model.rhs = b;
    cp.Model.lb = -Inf(size(f));
    cp.Model.ub = Inf(size(f));
    solve(cp)
%     fprintf('i %d\n',i);
% end