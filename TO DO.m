%% finish cbc upgrade CHECK LIFT AND PROJECT CUTS ON testMILP2

%% clp qp problem
% https://projects.coin-or.org/Clp/ticket/60
clc
prob = coinRead('crashqp.qps'); 


prob.H = tril(prob.H);
% prob.f = [1;1];


% prob.lb = [0;0];
% prob.ub = [1;1];

popts = [];
popts.display = 5;

clp(prob.H, prob.f, prob.A, prob.rl, prob.ru, prob.lb, prob.ub, popts)


%% cbc problems lift and project cuts
% https://projects.coin-or.org/Cbc/ticket/135


%% LEVMAR single inequality problem
% OPTI problems\ParameterEstimationProblem

%% MOSEK SDP support