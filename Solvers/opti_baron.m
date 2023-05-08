function [x,fval,exitflag,info] = opti_baron(nlprob,x0)
%OPTI_BARON Solve a LP/MILP/QP/QCQP/MIQP/MIQCQP/NLP/MINLP using BARON
%
%   [x,fval,exitflag,info] = opti_baron(nlprob,x0) solves optimization
%   problem specified by nlprob, together with optional initial guess x0.
%
%   THIS IS A WRAPPER FOR BARON

%   Copyright (C) 2013 Jonathan Currie (Control Engineering)

t = tic;
%Check required fields
if(~isfield(nlprob,'fun') || ~isfield(nlprob,'options'))
    error('You must use solve(optiObj) to generate a problem structure for this function');
end
%Check if we have a starting guess
if(nargin < 2 || isempty(x0))
    if(isfield(nlprob,'x0') && ~isempty(nlprob.x0))
        x0 = nlprob.x0;
    else
        x0 = [];
    end
end

% Run BARON
[x,fval,ef,inf,allsol] = baron(nlprob.fun,nlprob.A,nlprob.rl,nlprob.ru,nlprob.lb,nlprob.ub,...
                               nlprob.nlcon,nlprob.cl,nlprob.cu,nlprob.xtype,x0,nlprob.options);

%Collect Results
if(~isempty(inf))
    info.BBNodes = inf.BaR_Iterations;
end
info.Time = toc(t);
info.Algorithm = 'BARON: Branch And Reduce Optimization Navigator';
info.AllSol = allsol;
if(isempty(inf) || ~isfield(inf,'BARON_Status'))
    info.Status = 'Unknown';
    exitflag = -4;
else
    switch(ef)
        case 1
            exitflag = 1;
            info.Status = inf.Model_Status;
        case 2
            exitflag = -1;
            info.Status = inf.Model_Status;
        case 3
            exitflag = -2;
            info.Status = inf.Model_Status;
        case 4
            exitflag = -3;
            info.Status = inf.Model_Status;
        otherwise
            switch(lower(inf.BARON_Status))
                case {'max nodes in memory exceeded','max iterations exceeded','max cpu time exceeded'}
                    exitflag = 0;
                    info.Status = inf.BARON_Status;
                case 'run interrupted by user'
                    exitflag = -5;
                    info.Status = inf.BARON_Status;
                otherwise
                    exitflag = -4;
                    info.Status = inf.BARON_Status;                    
            end
    end
end