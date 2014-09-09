function [tstop,res] = testSolver(prob,solver,no,wbar)
%TESTSOLVER Perform a speed and accuracy test of an OPTI solver
%
%   [times,results] = testSolver(prob,solver) runs a series of tests 
%   across the selected solver given a problem type, prob, returning the 
%   execution times as well as a results structure.
%
%   [..] = test_Solver(prob,solver,no) runs 'no' number of problems
%   across each solver.

%   Copyright (C) 2011 Jonathan Currie (I2C2)

if(nargin < 2)
    error('You must supply the problem type and solver to testSolver');
end

%Setup default arguments
if(~exist('no','var'))
    no = [];
end
if(~exist('solver','var') || isempty(solver))
    solver = checkSolver(['best_' prob]);
end
if(~exist('wbar','var') || isempty(wbar))
    h = waitbar(0,['Solver: ' upper(solver)],'name','OPTI Solver Test','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    x0 = 0;
    n = 1;
    newH = 1;
else
    h = wbar.h;
    x0 = wbar.x0;
    n = wbar.n;
    newH = 0;
end
%Setup problem generation function
switch(lower(prob))
    case 'lp'
        testSet = @lp_prob;
    case 'milp'
        testSet = @milp_prob;
    case 'bilp'
        testSet = @bilp_prob;
    case 'qp'
        testSet = @qp_prob;
   case 'qcqp'
        testSet = @qcqp_prob;
    case 'miqp'
        testSet = @miqp_prob;
    case 'sdp'
        testSet = @sdp_prob;
    case 'nlp'
        testSet = @nlp_prob;
    case 'nls'
        testSet = @nls_prob;
    case 'minlp'
        testSet = @minlp_prob;
    otherwise
        delete(h);
        error('Unknown problem type: %s',prob);
end
%Setup no problems
if(isempty(no))
    no = testSet(); %default to maximum
elseif(no > testSet())
    warning('test:max','The maximum number of %s problems is %d',prob,testSet());
    no = testSet(); %cap at maximum
elseif(no < 1)
    no = 1;
end

%Preallocate
tstop = zeros(no,1);
res = struct('sol',{cell(no,1)},'acc',zeros(no,1),'fval',zeros(no,1),'ok',zeros(no,1),'msg',{cell(no,1)});

try
    %Initialize with two solves of problem 1
    prob = testSet(1);
    opts = optiset('solver',solver,'warnings','none');
    Opt = opti(prob,opts);
    solve(Opt); solve(Opt);
    %Now run tests
    for i = 1:no
        %Check Cancel
        if getappdata(h,'canceling')
            break;
        end
        [prob,~,fval] = testSet(i);
        Opt = opti(prob,opts);
        [x,sval,~,stat] = solve(Opt);
        res.sol{i} = x;
        res.acc(i) = abs(fval-sval);
        res.fval(i) = sval;
        try
            res.ok(i) = checkSol(Opt);
        catch ME
            if(strcmpi(ME.identifier,'MATLAB:innerdim') && res.acc(i) < 1e-6)
                res.ok(i) = true; %valid answer, but reduced the number of decision vars
            else
                res.ok(i) = false; 
            end
        end
        res.msg{i} = stat.Status;
        tstop(i) = stat.Time;
        %Update Waitbar
        waitbar(x0+i/no/n,h);    
    end
catch ME
    delete(h);
    if(exist('i','var'))
        error('Error testing problem %s - %s',num2str(i),ME.message);
    else
        error('Error initializing! - %s',ME.message);
    end
end

if(newH)
    delete(h);
end
