%% QCQP
clc
clear
prob = amplRead('testQCQP.mod',[],[],1);

funcs.objective = prob.fun;
funcs.gradient = prob.f;
funcs.constraints = prob.nlcon;
funcs.jacobian = prob.nljac;
funcs.jacobianstructure = prob.nljacstr;
funcs.hessian = prob.H;
funcs.hessianstructure = prob.Hstr;

opts.lb = prob.lb;
opts.ub = prob.ub;
opts.cl = prob.cl;
opts.cu = prob.cu;

opts.ipopt.print_level = 5;
% opts.ipopt.hessian_constant = 'yes'

[x,f] = ipopt(prob.x0,funcs,opts)


%%
clc
prob2 = amplRead('testQCQP.mod')
prob2.type = 'QCQP';
prob2.int.str = prob2.int;
prob2.sizes.nqc = length(prob2.qrl);

prob3 = QCQP2NLP(prob2,optiset)
prob3 = rowlin2nl(prob3,1,0)

funcs3.objective = prob3.fun;
funcs3.gradient = prob3.f;
funcs3.constraints = prob3.nlcon;
funcs3.jacobian = @(x) sparse(prob3.nljac(x));
funcs3.jacobianstructure = prob3.nljacstr;
funcs3.hessian = @(x,sigma,lambda) tril(prob3.H(x,sigma,lambda));
funcs3.hessianstructure = @() tril(prob3.Hstr());

opts = [];
opts.lb = prob3.lb;
opts.ub = prob3.ub;
opts.cl = prob3.cl;
opts.cu = prob3.cu;

opts.ipopt.print_level = 5;
% opts.ipopt.hessian_constant = 'yes'

[x,f] = ipopt(prob3.x0,funcs3,opts)

%%
clc
prob4 = amplRead('testQCQP.mod')
Opt = opti(prob4,optiset('solver','ipopt','display','iter'))

% Opt.nlprob.options.ipopt = [];
% Opt.nlprob.options.ipopt.print_level = 5;

[x,fval,ef,info] = solve(Opt)
