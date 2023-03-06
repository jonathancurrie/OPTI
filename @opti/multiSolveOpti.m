function [x,fval,exitflag,info] = multiSolveOpti(optObj,user_x0,ndivs,penalty,solveAll)
%MULTISOLVESOPTI an OPTI object
%
%   Called By opti multisolve

%   Very basic and naive implementation of a multi-start solver, using
%   pretty much an exhaustive search. However it has been very useful on
%   problems where an initial guess is hard to find!

%   Copyright (C) 2013 Jonathan Currie (Control Engineering)

%Number of 'best points' to keep during search
nbestpts = 10;
%Number of points to use for phase 2
nphase2pts = ceil(nbestpts/2);
%Number of divisions for phase 2 grid
ndivs2 = 2; %over written below if user specifies a vector
%Number of points to use for 'around x0' checking
nx0pts = 5;
%'around x0' radius factor
nx0rad = 0.05;

%If penalty is empty, set default
if(isempty(penalty)), penalty = 1e4; end    

%Get display level
verb = dispLevel(optObj.opts.display);
%Build a new OPTI object for use here, disabling output as required (note sense = 1 to prevent double inversion)
Opt = opti(optObj,'sense',1,'options',optiset(optObj.opts,'display','off','warnings','none','derivCheck','off'));

%Allocate input args
prob = Opt.prob; 
opts = Opt.opts;
%Initialization
where = '';
x_best = []; fval_best = inf; ef_best = []; info_best = []; i_best = 1;
titer = 0; tfeval = 0;

%Must have bounds!
if(isempty(prob.lb) || isempty(prob.ub) || any(isinf(prob.lb)) || any(isinf(prob.ub)))
    error('Multisolve only works for problems with finite rectangular bounds');
end

%If divs is empty, attempt to solve for a version that results in approx 100 sols, min 2 divisions
ndec = length(prob.lb);
if(isempty(ndivs))
    if(ndec==1)
        ndivs = 25;
    else
        ndivs = max(ceil(10^(2/ndec)),2);
    end
elseif(length(ndivs) == 2)
    ndivs2 = ndivs(2);
    if(ndivs2 < 2), nphase2pts = 0; end
    ndivs = ndivs(1);
end
%Determine number of points we are going to check
total_pts = ndivs^ndec;
%For search problems, limit nop based on total_pts
nbestpts = min(nbestpts,total_pts);

%Check division input
if(ndivs(1) == 1), error('You must have more than one division to use this function'); end
if(ndivs(1) > 1e4), error('The maximum number of divisions is 1e4'); end

%Check we have finite bounds
if(isempty(prob.lb) || any(isinf(prob.lb))), error('Multisolve currently requires finite lower bounds on all variables'); end
if(isempty(prob.ub) || any(isinf(prob.ub))), error('Multisolve currently requires finite upper bounds on all variables'); end

%Start Timer
t = tic;
%Build bound checking vectors
cvec = cell(ndec,1);
for i = 1:ndec
    cvec{i} = linspace(prob.lb(i),prob.ub(i),ndivs+1);
end

%Display header
if(verb)
    fprintf('------------------------------------------------------\n');
    fprintf('OPTI Multi-Start Solver using %s (beta)\n',upper(opts.solver));    
end

%Construct bound selection matrix for phase 1
lbnd = zeros(total_pts,ndec); lbnd2 = [];
ubnd = zeros(total_pts,ndec); ubnd2 = [];
col = ones(1,ndec+1);
for i = 1:total_pts
    for j = 1:ndec        
        lbnd(i,j) = cvec{j}(col(j));
        ubnd(i,j) = cvec{j}(col(j)+1);
    end
    col(1) = col(1) + 1;
    ind = col > ndivs;
    while(any(ind))
        idx = find(ind);
        col(idx) = 1;
        col(idx+1) = col(idx+1) + 1;
        ind = col > ndivs;
    end
end

if(verb)
    if(~solveAll)
        fprintf('\n Phase 1: Exhaustively Searching %d bounded regions (%d divisions per variable):\n',total_pts,ndivs);
    else
        fprintf('\n Exhaustively Solving %d bounded regions (%d divisions per variable):\n',total_pts,ndivs);
    end
end

%Begin searching bounded regions
x_search = cell(nbestpts,1); x_search2 = []; 
lb_search = cell(nbestpts,1); ub_search = cell(nbestpts,1); 
f_search = Inf(nbestpts,1); f_search2 = [];  
for i = 1:total_pts
    %Assign Bounds for this iter
    lb = lbnd(i,:)'; ub = ubnd(i,:)';
    %Make x0
    x0 = (ub-lb)./2 + lb;
    
    %Search Every Point
    if(~solveAll)    
        [x_search,f_search,lb_search,ub_search] = searchPoint(prob,x0,lb,ub,penalty,verb,i,total_pts,x_search,f_search,lb_search,ub_search);        
    %Else Solve Every Point
    else
        [x_best,fval_best,ef_best,info_best,i_best,where,titer,tfeval] = solvePoint(Opt,lb,ub,x0,verb,i,total_pts,titer,tfeval,x_best,fval_best,ef_best,info_best,i_best,where,'Exhaustive Solve');
    end   
end

%If we are not solving every problem, now solve all best points within their bounds, save best solution then enter phase 2
if(~solveAll)
    if(verb), fprintf('\n Solving Problem in the Best %d Phase 1 Regions:\n',nbestpts); end
    for i = 1:length(x_search)
        %Solve the problem within the search region
        [x_best,fval_best,ef_best,info_best,i_best,where,titer,tfeval] = solvePoint(Opt,lb_search{i},ub_search{i},x_search{i},verb,i,nbestpts,titer,tfeval,x_best,fval_best,ef_best,info_best,i_best,where,'Phase 1');     
    end  
    
    %Enter Phase 2, take noPhase2 best boxes and divide into divs2 divisions    
    noPhase2 = min(nphase2pts,length(lb_search));
    if(noPhase2)
        total_pts2 = ndivs2^ndec;
        if(verb)
            fprintf('\n Phase 2: Exhaustively Searching %d Sub-Regions Within the %d Best Phase 1 Regions using %d Divisions...\n',total_pts2*noPhase2,noPhase2,ndivs2);
        end
        cvec2 = cell(ndec,noPhase2);
        for i = 1:noPhase2
            for j = 1:ndec
                cvec2{j,i} = linspace(lb_search{i}(j),ub_search{i}(j),ndivs2+1);
            end
        end
        lbnd2 = zeros(total_pts2*noPhase2,ndec); 
        ubnd2 = zeros(total_pts2*noPhase2,ndec);     
        for p = 1:noPhase2
            col = ones(1,ndec+1);
            for i = 1:total_pts2
                for j = 1:ndec        
                    lbnd2(i + total_pts2*(p-1),j) = cvec2{j,p}(col(j));
                    ubnd2(i + total_pts2*(p-1),j) = cvec2{j,p}(col(j)+1);
                end
                col(1) = col(1) + 1;
                ind = col > ndivs2;
                while(any(ind))
                    idx = find(ind);
                    col(idx) = 1;
                    col(idx+1) = col(idx+1) + 1;
                    ind = col > ndivs2;
                end
            end
        end
        %Begin searching bounded regions
        x_search2 = cell(nbestpts,1); lb_search2 = cell(nbestpts,1); ub_search2 = cell(nbestpts,1); f_search2 = Inf(nbestpts,1);
        titer = 0; tfeval = 0; t = tic;
        for i = 1:total_pts2*noPhase2
            %Assign Bounds for this iter
            lb = lbnd2(i,:)'; ub = ubnd2(i,:)';
            %Make x0
            x0 = (ub-lb)./2 + lb;            
            %Search Point
            [x_search2,f_search2,lb_search2,ub_search2] = searchPoint(prob,x0,lb,ub,penalty,verb,i,total_pts2*noPhase2,x_search2,f_search2,lb_search2,ub_search2);  
        end

        if(verb), fprintf('\n Solving Problem in the Best %d Phase 2 Regions:\n',nbestpts); end
        for i = 1:length(x_search2)
            %Solve the problem within the search region
            [x_best,fval_best,ef_best,info_best,i_best,where,titer,tfeval] = solvePoint(Opt,lb_search2{i},ub_search2{i},x_search2{i},verb,i,length(x_search2),titer,tfeval,x_best,fval_best,ef_best,info_best,i_best,where,'Phase 2');     
        end 
        
        %Solve Problem within min and max bounds of best phase 2 divisions
        llb = Inf(ndec,1); gub = -Inf(ndec,1);
        for i = 1:length(lb_search2)
            for j = 1:ndec
                llb(j) = min([llb(j) lb_search2{i}(j) ub_search2{i}(j)]);
                gub(j) = max([gub(j) lb_search2{i}(j) ub_search2{i}(j)]);
            end
        end
        if(verb), fprintf('\n Solving Problem within Bounds of the Best Phase 2 Regions\n'); end
        [x_best,fval_best,ef_best,info_best,i_best,where,titer,tfeval] = solvePoint(Opt,llb,gub,(gub-llb)./2 + llb,verb,1,1,titer,tfeval,x_best,fval_best,ef_best,info_best,i_best,where,'Phase 2 Bounds');     
    end
end

%Solve from user supplied start point
if(~isempty(user_x0))
    if(verb), fprintf('\n Phase 3: Searching Problem around User Supplied Start Point:\n'); end
    %Create search cells
    x_searchx0 = cell(nx0pts+1,1);
    f_searchx0 = Inf(nx0pts+1,1); 
    %Variance of search points about user x0
    var = (prob.ub-prob.lb).*nx0rad;
    %Search user point
    [x_searchx0,f_searchx0] = searchPoint(prob,user_x0,[],[],penalty,verb,1,nx0pts+1,x_searchx0,f_searchx0); 
    %Search random points
    for i = 2:nx0pts+1
        %Create random search points within bounds
        x0 = randnbnd(user_x0,var,prob.lb,prob.ub);
        %Search Point
        [x_searchx0,f_searchx0] = searchPoint(prob,x0,[],[],penalty,verb,i,nx0pts+1,x_searchx0,f_searchx0); 
    end
    if(verb), fprintf('\n Solving Problem from the Best 3 Phase 3 Points:\n'); end
    %Solve Problem
    [x_best,fval_best,ef_best,info_best,i_best,where,titer,tfeval] = solvePoint(Opt,prob.lb,prob.ub,user_x0,verb,1,3,titer,tfeval,x_best,fval_best,ef_best,info_best,i_best,where,'User x0');
    for i = 1:2
        [x_best,fval_best,ef_best,info_best,i_best,where,titer,tfeval] = solvePoint(Opt,prob.lb,prob.ub,user_x0,verb,i+1,3,titer,tfeval,x_best,fval_best,ef_best,info_best,i_best,where,'Point Around User x0');
    end
else
    x_searchx0 = [];
    f_searchx0 = [];
end

%Check we have a solution so far
if(isempty(x_best))
    error('OPTI Multi-Solve did not find a feasible/successful solution.\n\nPlease increase the number of search points via the third argument to multisolve()%s','.');
end

%Check best solution in original problem
if(verb)
    if(~solveAll)
        fprintf('\nSolving Original Problem from the Best Point found...\n');
    else
        fprintf('\nRe-solving Original Problem from the Best Point found...\n');
    end
end
%Now Resolve from the best point, given original problem
[x,fval,exitflag,info] = solveOpti(optObj,x_best);
%Sum Iterations + FuncEvals
titer = incField(info,'Iterations',titer);
tfeval = incField(info,'FuncEvals',tfeval);
%Print Results
if(verb)
    status = getStatusMsg(exitflag);
    if(~isfield(info,'Iterations'))
        fprintf('Final Solution: Fval %12.5g, Status %s\n',fval,status);
    else
        fprintf('Final Solution: Fval %12.5g, Iter %4d, Status %s\n',fval,info.Iterations,status);
    end
end
%Check if final solution is actually better than best solution
if(fval_best < fval)
    x = x_best;
    fval = fval_best;
    exitflag = ef_best;
    info = info_best;
    if(verb)
        fprintf('\nFinal Solution is worse than best found, returning best found (%s, Run %d)\n',where,i_best);
    end
end
    
%Modify Info Status
info.Time = toc(t);
if(titer), info.Iterations = titer; end
if(tfeval), info.FuncEvals = tfeval; end

%Save Search Area
optObj.prob.multi.lbnd = lbnd;
optObj.prob.multi.ubnd = ubnd;
optObj.prob.multi.lbnd2 = lbnd2;
optObj.prob.multi.ubnd2 = ubnd2;
optObj.prob.multi.x_search = x_search;
optObj.prob.multi.f_search = f_search;
optObj.prob.multi.x_search2 = x_search2;
optObj.prob.multi.f_search2 = f_search2;
optObj.prob.multi.x_searchx0 = x_searchx0;
optObj.prob.multi.f_searchx0 = f_searchx0;

if(verb)
    fprintf('------------------------------------------------------\n');
end



%Search a Point, print as required
function [x_search,f_search,lb_search,ub_search] = searchPoint(prob,x0,lb,ub,penalty,verb,i,total_pts,x_search,f_search,lb_search,ub_search)        
if(nargin < 12), ub_search = []; end
if(nargin < 11), lb_search = []; end
%Evaluate point
[fpen,f,c] = evalPoint(prob,x0,penalty);
%Decide whether we save this point
idx = fpen < f_search; 
bstr = 1;
if(any(idx))
    %Find first index we are better than, insert, and shift remainder
    fidx = find(idx); fidx = fidx(1);
    %Insert
    f_search = [f_search(1:fidx-1);fpen;f_search(fidx:end-1)];
    x_search = [x_search(1:fidx-1);x0;x_search(fidx:end-1)];
    if(~isempty(lb_search))
        lb_search = [lb_search(1:fidx-1);lb;lb_search(fidx:end-1)];
        ub_search = [ub_search(1:fidx-1);ub;ub_search(fidx:end-1)];
    end
    %Indicate new solution found
    bstr = 2;
end
if(verb)            
    if(~isempty(c))
        fprintf(bstr,'Point %3d of %3d: FPenalty %12.5g, Fval %12.5g, ConViol %12.5g\n',i,total_pts,fpen,f,c);
    else
        fprintf(bstr,'Point %3d of %3d: Fval %12.5g\n',i,total_pts,f);
    end
end

%Solve a Point, print as required
function [x_best,fval_best,ef_best,info_best,i_best,where,titer,tfeval] = solvePoint(Opt,lb,ub,x0,verb,i,total_pts,titer,tfeval,x_best,fval_best,ef_best,info_best,i_best,where,phase)       
%Save new bounds into OPTI object
Opt.prob.lb = lb; Opt.prob.ub = ub;
Opt.nlprob.lb = lb; Opt.nlprob.ub = ub;
Opt.nlprob.options.lb = lb; Opt.nlprob.options.ub = ub;
bstr = 1;
%Solve Problem
try
    %Solve Problem
    [xi,fi,ei,ii] = solveOpti(Opt,x0); 
    %Save results into object
    Opt.sol = xi; Opt.obj = fi; Opt.ef = ei; Opt.info = ii;
    %If best, save it
    if(fi < fval_best)
        %Check point is actually feasible
        if(checkOptiSol(Opt,1e-4)) %slightly relaxed tol
            x_best = xi; fval_best = fi; ef_best = ei; info_best = ii; i_best = i; bstr = 2; where = phase;
        end
    end        
    %Sum Iterations + FuncEvals
    titer = incField(ii,'Iterations',titer);
    tfeval = incField(ii,'FuncEvals',tfeval);
catch ME
    fi = NaN; ei = NaN; ii = NaN;
    optiwarn('OPTI:SolverError',ME.message);
end
%Display Iteration if requested
if(verb)
    status = getStatusMsg(ei);
    if(~isfield(ii,'Iterations'))
        fprintf(bstr,'Run %3d of %3d: Fval %12.5g, Status %s\n',i,total_pts,fi,status);
    else
        fprintf(bstr,'Run %3d of %3d: Fval %12.5g, Iter %4d, Status %s\n',i,total_pts,fi,ii.Iterations,status);
    end
end
        
        

function [fpen,f,c] = evalPoint(prob,x0,penalty)
%Evaluate point
f = prob.objective(x0);
if(~isempty(prob.constraints))
    c = sum(prob.constraints(x0));    
    fpen = f + c*penalty;
else
    c = [];
    fpen = f;
end

function status = getStatusMsg(ef)
switch(ef)
    case 1
        status = 'OK';
    case 0
        status = 'Exceeded It/Time';
    case -1
        status = 'Infeasible';
    otherwise
        status = 'Error';
end

function val = incField(strct,field,val)
if(isfield(strct,field) && isnumeric(strct.(field)))
    val = val + strct.(field);
end

function x0 = randnbnd(mean,var,lb,ub)
ndec = length(lb);
%Clamp x0 outside bounds
il = mean < lb; mean(il) = lb(il);
iu = mean > ub; mean(iu) = ub(iu);
%Check if avg sitting on a bound, if so try move by var in
for i = 1:ndec
    if(mean(i) == lb(i))
        mean(i) = min(mean(i) + var(i),ub(i));
    elseif(mean(i) == ub(i))
        mean(i) = max(mean(i) - var(i),lb(i));
    end
end
%Now generate random numbers until all fall within bounds, or 100 iter
idx = true(ndec,1); iter = 1; x0 = zeros(ndec,1);
while(any(idx) && iter < 100)
    xnew = mean(idx) + var(idx).*randn(sum(idx),1);
    x0(idx) = xnew;
    idx = x0 < lb | x0 > ub;
    iter = iter + 1;
end
if(iter==100)
    %re clamp
    il = x0 < lb; x0(il) = lb(il);
    iu = x0 > ub; x0(iu) = ub(iu);
end


