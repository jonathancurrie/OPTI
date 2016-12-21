function [tstop,res] = testProblem(optObj,solvers)
%TESTPROBLEM  Test an OPTI Problem against all available solvers
%
%   [times,results] = testProblem(optObj) runs the supplied OPTI problem
%   against all available solvers for the given problem type.
%
%   [times,results] = testProblems(optObj,solvers) runs only the selected
%   solvers across the supplied problem.

%   Copyright (C) 2011 Jonathan Currie (IPL)


%Setup default arguments
if(~exist('solvers','var') || isempty(solvers))
    solvers = optiSolver(['all_' optObj.prob.type]);
end
if(~iscell(solvers))
    solvers = {solvers};
end

%Remove troublesome solvers
ind = strcmpi(solvers,'lipsol');
if(any(ind))
    solvers(ind) = [];
end

%Setup problem generation function (for warm up)
switch(lower(optObj.prob.type))
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
    case 'sdp'
        testSet = @sdp_prob;
    case 'miqp'
        testSet = @miqp_prob;
    case 'nlp'
        testSet = @nlp_prob;
        if(isempty(optObj.prob.x0))
            error('Please include x0 in optiprob when using this function on NLPs');
        end
    case 'minlp'
        testSet = @minlp_prob;
        if(isempty(optObj.prob.x0))
            error('Please include x0 in optiprob when using this function on MINLPs');
        end
    otherwise
        error('Unknown problem type: %s',optObj.prob.type);
end

%Remove troublesome solvers
if(strcmpi(optObj.prob.type,'lp'))
    ind = strcmpi(solvers,'SEDUMI'); solvers(ind) = []; %no good with neq = nvar
end

%Setup no runs
len = length(solvers);
%confirm solver exists!
for i = 1:len
    rub = optiSolver(solvers{i}); %#ok<NASGU>
end

%Setup wait bar
h = waitbar(0,['Solver: ' upper(solvers{1})],'name','OPTI Problem Test','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

%Preallocate
tstop = zeros(len,1);
res = struct('sol',{cell(len,1)},'fval',zeros(len,1),'ok',zeros(len,1),'msg',{cell(len,1)});

%Clear display warnings etc
opts = optiset(optObj.prob,'display','off','warnings','off');

try
    for i = 1:len
        %Check Cancel
        if getappdata(h,'canceling')
            break;
        end
        %Assign the new solver
        solver = solvers{i};
        opts = optiset(opts,'solver',solver);
        %Initialize with two solves of problem 1
        prob = testSet(1);
        Opt = opti(prob,opts);
        solve(Opt); solve(Opt);
        %Now run our problem
        Opt = opti(getProb(optObj),opts);
        [x,sval,~,stat] = solve(Opt);
        res.sol{i} = x;
        res.fval(i) = sval;
        res.ok(i) = checkSol(Opt);
        res.msg{i} = stat.Status;
        tstop(i) = stat.Time;
        %Update Waitbar
        waitbar(i/len,h,['Solver: ' upper(solvers{min(len,i+1)})]);    
    end
catch ME
    delete(h);
    rethrow(ME);
end

delete(h);


%Draw Table
fprintf('\nTest Problem [%s] Result\n',Opt.prob.type);
fprintf('--------------------------------------------------\n');
fprintf('  Solver    Status        Fval         Time\n');
for i = 1:len
    fprintf('%8s    %4s    %+12.4f     %1.4fs\n',solvers{i},status(res.ok(i)),res.fval(i),tstop(i));
end
fprintf('--------------------------------------------------\n');


function msg = status(ok)

if(ok)
    msg = 'OK';
else
    msg = 'FAIL';
end
