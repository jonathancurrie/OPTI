%% UNO Rosenbrock
clc
%Objective
obj = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
%Setup Options
opts = optiset('solver','ipopt','display','iter');
%Build & Solve
Opt = opti('obj',obj,'ndec',2,'options',opts)
x0 = [0 0]';
[x,fval,exitflag,info]= solve(Opt,x0)
%Plot
% plot(Opt,[],1)

%% NLP1 Hock & Schittkowski #71
clc
%Objective & Gradient
obj = @(x) x(1)*x(4)*sum(x(1:3)) + x(3);
grad = @(x) [ x(1)*x(4) + x(4)*sum(x(1:3));
              x(1)*x(4);
              x(1)*x(4) + 1;
              x(1)*sum(x(1:3)) ];          
%Linear Constraints
lb = ones(4,1);
ub = 5*ones(4,1);
%Nonlinear Constraints
nlcon = @(x) [ prod(x);
               sum(x.^2)];
nljac = @(x) [ prod(x)./x';
                2*x' ];          
nlrhs = [25 40]';
nle = [1 0]'; % (>=, ==)
%Setup Options
opts = optiset('solver','ipopt','warnings','on','display','iter','solverOpts',ipoptset('linear_solver','pardiso'));
%Build & Solve
Opt = opti('obj',obj,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'options',opts)
x0 = [1 5 5 1]';
[x,fval,exitflag,info]= solve(Opt,x0)

info.Lambda

%%
clc
prob = amplRead('opf_662bus'); 
opts = optiset('solver','ipopt','solverOpts',ipoptset('linear_solver','pardiso'),'display','off');
opts2 = optiset('solver','ipopt','solverOpts',ipoptset('linear_solver','ma57'),'display','off');
opts3 = optiset('solver','ipopt','solverOpts',ipoptset('linear_solver','mumps'),'display','off');
O = opti(prob,opts)
O2 = opti(prob,opts2);
O3 = opti(prob,opts3);

[x,f,e,i] = solve(O);
[x2,f,e,i2] = solve(O2);
[x3,f,e,i3] = solve(O3);

norm(x-x2)

i
i2
i3



%%
clc
clear funcs opts

x0 = [0; 0];

funcs.objective = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
funcs.gradient = @(x)[2*x(1)-400*x(1)*(x(2)-x(1)^2)-2,200*x(2)-200*x(1)^2];

opts.ipopt.print_level = 5;
opts.ipopt.hessian_approximation = 'limited-memory';
opts.ipopt.ma57_pivot_order = 2;
opts.ipopt.linear_solver = 'pardiso';

[x,f] = ipopt(x0,funcs,opts)

f.eval


%%

prob = amplRead('ch3.nl')
opts = optiset('solver','ipopt','display','iter','solverOpts',ipoptset('linear_solver','pardiso'));

Opt = opti(prob,opts)

x = solve(Opt)

asl('close')