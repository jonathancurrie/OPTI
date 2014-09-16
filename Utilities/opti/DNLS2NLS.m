function prob = DNLS2NLS(prob,opts)
%DNLS2NLS  Converts a Dynamic NLS problem into a NLS problem
% prob= DNLS2NLS(prob,opts)

%   Copyright (C) 2013 Jonathan Currie (I2C2)

if(~isfield(prob,'type') || isempty(prob.type))
    error('This function is not for user use, and should only be called from OPTI');
end

%Get warning level
warn = optiWarnLevel(opts.warnings);

%Get Derivative Checker Level
derCheck = false;
if(strcmpi(opts.derivCheck,'on')), derCheck = true; end

%Transpose as required
if(size(prob.odez0,2) > 1), prob.odez0 = prob.odez0'; end

%Get Options & Set Sizes
dopts = optidynset(opts.dynamicOpts);
dopts.n = length(prob.odez0);           %number of states
dopts.inz0 = isnan(prob.odez0);         %index of initial states to estimate
dopts.nz0 = sum(dopts.inz0);            %number of initial states to estimate
dopts.np = length(prob.x0)-dopts.nz0;   %number of parameters
dopts.inMeasIdx = 0; %logical index
%Set Maximum Integrator Time if no output function specified
if(~isfield(dopts.odeOpts,'OutputFcn') || isempty(dopts.odeOpts.OutputFcn))
    dopts.odeOpts = odeset(dopts.odeOpts,'OutputFcn',@(t,y,flag) odeMaxTime(t,y,flag,dopts.odeMaxTime));
else
    if(warn > 1), optiwarn('OPTI:NoOutFcn','ODE Output Function has already been specified by the user, ODE max time cannot be set'); end
end
%Get Initial Time 
t0 = dopts.initialT;

%If user supplied a stateIndex, we may not be comparing all states with measured data
if(~isempty(dopts.stateIndex))
    if(isnumeric(dopts.stateIndex))
        if(any(dopts.stateIndex > dopts.n) || any(dopts.stateIndex < 1))
            error('State Index vector contains indices < 1 or > number of states');
        end
        if(length(dopts.stateIndex) ~= length(unique(dopts.stateIndex)))
            error('There appear to be repeated entries in the State Index vector');
        end
        if(any(sort(dopts.stateIndex) ~= dopts.stateIndex))
            error('The State Index vector must be sorted');
        end
        dopts.inz = false(1,dopts.n);
        dopts.inz(dopts.stateIndex) = true;        
    elseif(islogical(dopts.stateIndex))
        if(length(dopts.stateIndex) ~= dopts.n)
            error('When specifying a logical state index vector it must have the same number of elements as states');
        end
        dopts.inz = dopts.stateIndex;
    end
else
    dopts.inz = true(1,dopts.n);
end
G = false(dopts.n,dopts.np+dopts.nz0);
G(dopts.inz,:) = true;
dopts.inJacz = [false(1,dopts.n) reshape(G==true,1,numel(G))];

%HACK if user supplied initialT then we must use the cell based code
if(~isempty(t0) && ~iscell(prob.xdata) && ~iscell(prob.ydata))
    %Reshape ydata into matrix 
    ydata = reshape(prob.ydata,length(prob.xdata),sum(dopts.inz));    
    xd = cell(sum(dopts.inz),1); yd = xd;
    for i = 1:sum(dopts.inz)
        xd{i} = prob.xdata; %just copies
        yd{i} = ydata(:,i); %extract each measurement vector into a cell
    end
    prob.xdata = xd;
    prob.ydata = yd;
end

%Check ODE Function given above sizes
try
    prob.ode(1,ones(dopts.n,1),ones(dopts.np,1));
catch ME
    throwAsCaller(MException('OPTI:ODE_ERROR','There was an error evaluating your ODE function. This is normally due to an incorrectly specified x0 and/or z0.\n\nCheck the following error for more information:\n\n%s',ME.message));
end

%Check for xdata/ydata cell arrays (indicating multiple sampling rates between states)
if(xor(iscell(prob.xdata),iscell(prob.ydata)))
    throwAsCaller(MException('OPTI:ODE_ERROR','Both xdata and ydata must be a cell array if one is specified as a cell'));
end
if(iscell(prob.xdata))
    %Check lengths
    if(length(prob.xdata) ~= length(prob.ydata))
        error('Cell arrays xdata and ydata should contain the same number of cells');
    end
    %Check number of cells in xdata/ydata = no states (one cell per ode)
    if(length(prob.xdata) ~= sum(dopts.inz))
        error('When solving parameter estimation problems with multiple sampling rates, each data cell array should contain measurements for each ODE (i.e. xdata = ydata = cell(nstates,1))');
    end
    %If starting times are not all equal AND we are not estimating the
    %initial states of the LAST starting times, initial conditions the user
    %supplied MIGHT be incorrect
    startT = zeros(length(prob.xdata),1);
    for i = 1:length(startT)
        startT(i) = min(prob.xdata{i}); %use min incase t is not in order
    end
    %Augment artifical start time   
    if(length(unique(startT)) > 1)
        minT = min([t0;startT]);
        %Check if we estimating initial conditions if startTs > minT
        probInd = startT > minT & ~dopts.inz0(dopts.inz);
        stateNos = 1:dopts.n; stateNos = stateNos(dopts.inz);
        if(any(probInd))
            %Generate warning string
            str = sprintf('OPTI has detected a possible formulation error in the parameter estimation problem:\n\n');
            str = sprintf('%sThe following states have measurements starting later than t = %g:\n',str,minT);
            ind = find(probInd);
            for i = 1:length(ind)
                str = sprintf('%s State %d [first point at t = %g, z0[t=%g] = %g?]\n',str,stateNos(ind(i)),startT(ind(i)),minT,prob.odez0(stateNos(ind(i))));
            end
            str = sprintf('%s\nDue to one or more states having measurements starting at t = %g (or Initial Time set), OPTI will use',str,minT);
            str = sprintf('%s this as the ODE integration start time. This means your initial conditions (z0) for the above',str);
            str = sprintf('%s states must also be specified at t = %g. If this is not the case, please change these states''',str,minT);
            str = sprintf('%s initial conditions to being estimated as well (set corresponding z0 elements to NaN), to avoid this warning and possible solution errors.\n',str);
            if(warn), optiwarn('OPTI:DNLS_MultiRateIC',str); end
        end
    end        
    
    %Build Unique Time Vector, checking transpose and lengths as we go
    for i = 1:length(prob.xdata)
        if(length(prob.xdata{i}) ~= length(prob.ydata{i}))
            error('Fitting data in cell %d does not have a time stamp associated with each data point (lengths are not equal)',i);
        end
        if(~isempty(prob.weighting) && length(prob.weighting{i}) ~= length(prob.ydata{i}))
            error('Fitting data in cell %d does not have a weight associated with each data point (lengths are not equal)',i);
        end
        if(size(prob.xdata{i},1) > size(prob.xdata{i},2))
            prob.xdata{i} = prob.xdata{i}';
        end
    end
    ts = [prob.xdata{:}];
    if(size(ts,2) > 1), ts = ts'; end
    try
        tunq = unique([t0;ts],'legacy');
    catch
        tunq = unique([t0;ts]);
    end        
    
    %Determine if we have any repeated data points, and how many in each state (tedious, but required)
    noRepStates = zeros(1,sum(dopts.inz));
    for i = 1:length(tunq)
        for j = 1:length(prob.xdata)
            s = sum(prob.xdata{j} == tunq(i));
            if(s > 1)
                noRepStates(j) = max(noRepStates(j),s)-1;
            end
        end
    end
    
    %If any repeated points, bit of extra work to do here
    if(sum(noRepStates))
        %Create Fitting Data Matrix To Determine Indices
        fdata = NaN(length(tunq),length(prob.ydata)+sum(noRepStates));    
        wdata = NaN(length(tunq),length(prob.ydata)+sum(noRepStates)); wts = ~isempty(prob.weighting);
        for i = 1:length(tunq)
            for j = 1:length(prob.xdata)
                ind = prob.xdata{j} == tunq(i);
                if(any(ind))
                    %Determine which column(s) to insert data into
                    if(j > 1)
                        start = j + sum(noRepStates(1:j-1));                        
                    else
                        start = 1;
                    end
                    stop = start + sum(ind)-1;  %may not have all repeated measurements at the same time.... %noRepStates(j);
                    %Insert our data
                    fdata(i,start:stop) = prob.ydata{j}(ind);
                    if(wts)
                        wdata(i,start:stop) = prob.weighting{j}(ind);
                    end
                end
            end
        end    
        
        %Create index matrix to determine which states (and repeated measurement states) are in which columns        
        idxState = zeros(1,size(fdata,2));
        col = 1;
        for i = 1:length(noRepStates)
            for j = 0:noRepStates(i)
                idxState(col) = i;
                col = col + 1;
            end
        end
        
        %We now know where our measurements sit, now we must create an indexing vector using doubles
        nfdata = ~isnan(fdata);
        idx = find(nfdata(:,1)); len = size(nfdata,1);
        for i = 2:size(fdata,2)
            colidx = find(nfdata(:,i))+len*(idxState(i)-1);
            idx = [idx;colidx]; %#ok<AGROW> %should be fine in terms of speed, not many columns expected and large chunks
        end
        
        %Save Existing Values
        prob.misc.xdata_orig = prob.xdata;
        prob.misc.ydata_orig = prob.ydata;
        prob.misc.weighting_orig = prob.weighting;
        %Save new indexed xdata and ydata
        prob.xdata = tunq(:);
        ydata = fdata(~isnan(fdata));
        prob.ydata = ydata(:);
        if(wts)
            weights = wdata(~isnan(wdata));
            prob.weighting = weights(:);
        end
        %Remember our index measurement logical array will have to be doubles as we must repeat some outputs in our (yhat - ydata) calculation
        dopts.inMeas = idx; dopts.inMeasIdx = 2; %element based
        dopts.nmeas = length(dopts.inMeas);
        dopts.selectedMeas = true;
        %Build Jacobian Index Matrix
        dopts.inJacMeas = zeros(length(dopts.inMeas),dopts.np+dopts.nz0);
        dopts.inJacMeas(:,1) = dopts.inMeas;
        for i = 2:dopts.np+dopts.nz0
            dopts.inJacMeas(:,i) = dopts.inMeas + max(dopts.inMeas)*(i-1);
        end

    %Standard multiple-rate problem
    else    
        %Create Fitting Data Matrix To Determine Indices
        fdata = NaN(length(tunq),length(prob.ydata));
        wdata = NaN(length(tunq),length(prob.ydata)); wts = ~isempty(prob.weighting);
        for i = 1:length(tunq)
            for j = 1:length(prob.xdata)
                ind = prob.xdata{j} == tunq(i);
                if(any(ind))
                    fdata(i,j) = prob.ydata{j}(ind);
                    if(wts)
                        wdata(i,j) = prob.weighting{j}(ind);
                    end
                end
            end
        end
        dopts.inMeas = ~isnan(fdata);
        dopts.inJacMeas = repmat(dopts.inMeas,1,dopts.np+dopts.nz0);
        dopts.selectedMeas = true;
        %Save Old xdata, ydata
        prob.misc.xdata_orig = prob.xdata;
        prob.misc.ydata_orig = prob.ydata;
        prob.misc.weighting_orig = prob.weighting;
        %Save new indexed xdata and ydata
        prob.xdata = tunq(:);
        prob.ydata = fdata(dopts.inMeas);
        if(wts), prob.weighting = wdata(dopts.inMeas); end
        %Determine number of measurements (in all states)
        dopts.nmeas = sum(sum(dopts.inMeas));
    end
%Normal problem, every measurement at the same sampling intervals    
else   
    %Check User Input
    if(length(prob.xdata)*length(prob.odez0(dopts.inz)) ~= length(prob.ydata))
        error(['When solving Dynamic Parameter Estimation problems you must have the same number of points in '...
               'ydata as length(xdata)*nstates (or length(xdata) * the number of indexed states).\n\nlength(ydata) = %d [expected %d]\n'],...
               length(prob.ydata),length(prob.xdata)*length(prob.odez0(dopts.inz)));
    end
    
    %Order time stamps
    [prob.xdata,idx] = sort(prob.xdata);
    %Reshape ydata into matrix 
    ydata_old = reshape(prob.ydata,length(prob.xdata),sum(dopts.inz));
    %Order ydata rows based on time
    ydata_old = ydata_old(idx,:);
    %Save New Sorted Ydata    ;
    prob.ydata = ydata_old(:);
    %If we have a weighting matrix, reorder as above
    if(~isempty(prob.weighting))
        weights_old = reshape(prob.weighting,length(prob.xdata),sum(dopts.inz));
        weights_old = weights_old(idx,:);
        prob.weighting = weights_old(:);
    end        
    
    %Check for repeated measurements
    try
        [tunq,~,ic] = unique(prob.xdata,'legacy');
    catch
        [tunq,~,ic] = unique(prob.xdata);
    end
    %If lengths are different, we have repeated points
    if(length(tunq) ~= length(prob.xdata))
        %Save Existing Values
        prob.misc.xdata_orig = prob.xdata;
        prob.misc.ydata_orig = ydata_old;
        prob.misc.weighting_orig = prob.weighting;
        %Save unique values of time (otherwise ode int gets upset)
        prob.xdata = tunq;
        %Our index measurement logical array will now have to be doubles as we must repeat some outputs in our (yhat - ydata) calculation
        dopts.inMeas = ic; dopts.inMeasIdx = 1; %row based
        dopts.nmeas = length(dopts.inMeas)*sum(dopts.inz);
        dopts.selectedMeas = true;
        %Build Jacobian Index Matrix
        dopts.inJacMeas = zeros(length(dopts.inMeas),(dopts.np+dopts.nz0)*sum(dopts.inz));
        dopts.inJacMeas(:,1) = dopts.inMeas;
        for i = 2:(dopts.np+dopts.nz0)*sum(dopts.inz)
            dopts.inJacMeas(:,i) = dopts.inMeas + max(dopts.inMeas)*(i-1);
        end
    else
        %Normal Logical Index
        dopts.inMeas = true(length(prob.ydata),1);
        %Determine number of measurements (in all states)
        dopts.nmeas = sum(sum(dopts.inMeas));
        dopts.selectedMeas = false;
    end    
end

%Must have at least 3 points in t to avoid ode integrator interpreting as tspan
if(length(prob.xdata) < 3)
    error('You must supply at least 3 measurement points to solve a DNLS');
end

%Save number of time points
dopts.nT = length(prob.xdata);

%Decide what sensitivity strategy to use
if(isempty(dopts.sensitivity))
    if(~isempty(dopts.dfdz) || ~isempty(dopts.dfdp))
        dopts.sensitivity = 'User';
    else
        dopts.sensitivity = 'ND';
    end
end

%Check User Derivatives (if supplied) and/or Setup Derivative Estimation Method if required
z0 = prob.odez0; x0 = prob.x0;
if(dopts.nz0) %replace NaNs with initial guesses
    z0(dopts.inz0) = prob.x0(dopts.np+1:end);
    x0 = prob.x0(1:dopts.np);
end
if(strcmpi(dopts.sensitivity,'user'))
    if(~isempty(dopts.dfdz))
        if(nargin(dopts.dfdz) ~= 3), error('ODE dfdz must be a function which accepts 3 arguments (t,z,p)'); end  
        if(derCheck)
            z = @(z) prob.ode(1,z,x0);
            dfdz = @(z) dopts.dfdz(1,z,x0);
            optiDerCheck(z,dfdz,z0,'ODE DFDZ',warn);
        end
        dopts.dfdz_method = 2;
    else
        if(warn>1), optiwarn('OPTI:NoUserDer','You have specified to use user supplied derivatives for DFDZ, but the corresponding function is empty.\nUsing ND instead.'); end
        dopts.dfdz_method = 0;
    end
else
    %Setup derivative estimation method 
    if(strcmpi(dopts.sensitivity,'cs'))
        dopts.dfdz_method = 3;
    elseif(strcmpi(dopts.sensitivity,'ad'))
    	dopts.dfdz_method = 1;
    else
        dopts.dfdz_method = 0;
    end
end
if(strcmpi(dopts.sensitivity,'user'))
    if(~isempty(dopts.dfdp))
        if(nargin(dopts.dfdp) ~= 3), error('ODE dfdp must be a function which accepts 3 arguments (t,z,p)'); end
        if(derCheck)
            p = @(p) prob.ode(1,z0,p);
            dfdp = @(p) dopts.dfdp(1,z0,p);
            optiDerCheck(p,dfdp,x0,'ODE DFDP',warn);
        end
        dopts.dfdp_method = 2;
    else
        if(warn>1), optiwarn('OPTI:NoUserDer','You have specified to use user supplied derivatives for DFDP, but the corresponding function is empty.\nUsing ND instead.'); end
        dopts.dfdp_method = 0;
    end
else
    %Setup derivative estimation method
    if(strcmpi(dopts.sensitivity,'cs'))
        dopts.dfdp_method = 3;
    elseif(strcmpi(dopts.sensitivity,'ad'))
        dopts.dfdp_method = 1;
    else
        dopts.dfdp_method = 0;
    end
end


%Assign Objective
dopts.reqGrad = false;
prob.fun = @(theta,tm) odeEstim(prob.ode,prob.odez0,tm,theta,dopts);
prob.misc.fitFun = prob.fun;

%Assign Gradient (if required)
if(~strcmpi(dopts.sensitivity,'none'))
    %Setup Initial Sensitivity
    Sp0 = zeros(dopts.n,dopts.np);
    if(dopts.nz0)
        Sz0 = eye(dopts.n,dopts.n); 
        Sz0(:,~dopts.inz0) = [];
    else
        Sz0 = [];
    end
    dopts.S0 = [Sp0 Sz0];
    
    %Add Derivative Calculation
    dopts.reqGrad = true;
    prob.f = @(theta) odeEstim(prob.ode,prob.odez0,prob.xdata,theta,dopts);
    prob.misc.fitGrad = prob.f;
end

function f = odeEstim(odeFun,z0,tm,theta,dopts)
%Optimizer Callback Function   

%Replace z0 with optimizer guess if required
if(dopts.nz0)
    z0(dopts.inz0) = theta(dopts.np+1:end);
end

%Check if we need gradient
if(dopts.reqGrad)
    ode = @(t,x) odeSens(t,x,odeFun,theta(1:dopts.np),dopts);        
    %Augment Sensitivity initial states
    Z0 = [z0;dopts.S0(:)];
else
    ode = @(t,x) odeFun(t,x,theta(1:dopts.np)); %original function
    Z0 = z0;   
end

%Use selected integrator
switch(dopts.integrator)
    case 'ode45'
        [~,f] = ode45(ode,tm,Z0,dopts.odeOpts);        
    case 'ode15s'
        [~,f] = ode15s(ode,tm,Z0,dopts.odeOpts);   
    case 'ode23s'
        [~,f] = ode23s(ode,tm,Z0,dopts.odeOpts); 
    case 'ode23'
        [~,f] = ode23(ode,tm,Z0,dopts.odeOpts);
    case 'ode23t'
        [~,f] = ode23t(ode,tm,Z0,dopts.odeOpts); 
    case 'ode23tb'
        [~,f] = ode23tb(ode,tm,Z0,dopts.odeOpts); 
    case 'ode15i'
        [~,f] = ode15i(ode,tm,Z0,dopts.odeOpts); 
    case 'ode113'
        [~,f] = ode113(ode,tm,Z0,dopts.odeOpts); 
    otherwise
        error('Integrator ''%s'' not implemented yet',dopts.integrator);
end

if(dopts.reqGrad) 
    %Check we got the correct sized vector, otherwise integrator may have failed
    [r,c] = size(f);
    if(r~=dopts.nT || c~=dopts.n*(1+dopts.np+dopts.nz0))
        fprintf(2,'  INTEGRATOR ERROR IN GRADIENT EVALUATION - OPTI IS RETURNING 1e6 AT THIS POINT\n');
        f = 1e6*ones(dopts.nT,dopts.n*(1+dopts.np+dopts.nz0));
    end
    %Index Sensitivity States (n+1) & States we are comparing to measured data
    f = f(:,dopts.inJacz);
    %Index Measurements within the Jacobian we will compare to
    if(dopts.selectedMeas)
        f = f(dopts.inJacMeas);
    end
    %Reshape into standard Jacobian
    f = reshape(f,dopts.nmeas,dopts.nz0+dopts.np);    
else   
    %Check we got the correct sized vector, otherwise integrator may have failed
    [r,c] = size(f);
    if(r~=dopts.nT || c~=dopts.n)
        fprintf(2,'  INTEGRATOR ERROR IN OBJECTIVE EVALUATION - OPTI IS RETURNING 1e6 AT THIS POINT\n');
        f = 1e6*ones(dopts.nT,dopts.n);
    end
    %Index States we are comparing to measured data
    f = f(:,dopts.inz);
    %Index Measurements within those states we will compare to
    switch(dopts.inMeasIdx)
        case 0 %logical
            f = f(dopts.inMeas);
        case 1 %row
            f = f(dopts.inMeas,:);
            f = f(:);
        case 2 %element 
            f = f(dopts.inMeas);
    end
end



function zdot = odeSens(t,z,ode,theta,dopts)
%ODE with Integrated Sensitivity Function 
%Get dfdz
switch(dopts.dfdz_method)
    case 0 %numerical differentiation
        dfdz = mklJac(@(z) ode(t,z,theta),z(1:dopts.n));
    case 1 %automatic differentiation
        dfdz = autoJac(@(z) ode(t,z,theta),z(1:dopts.n));
    case 2 %user supplied
        dfdz = dopts.dfdz(t,z,theta);
    case 3 %complex step
        dfdz = cstepJac(@(z) ode(t,z,theta),z(1:dopts.n));
end
%Get dfdp
if(dopts.np > 0)
    switch(dopts.dfdp_method)
        case 0 %numerical differentiation
            dfdp = mklJac(@(theta) ode(t,z(1:dopts.n),theta),theta);
        case 1 %automatic differentiation
            dfdp = autoJac(@(theta) ode(t,z(1:dopts.n),theta),theta);
        case 2 %user supplied
            dfdp = dopts.dfdp(t,z,theta);
        case 3 %complex step
            dfdp = cstepJac(@(theta) ode(t,z(1:dopts.n),theta),theta);
    end
else
   dfdp = []; 
end
%Complete sensitivity differential equation
S = dfdz*reshape(z(dopts.n+1:end),dopts.n,dopts.nz0+dopts.np) + [dfdp zeros(dopts.n,dopts.nz0)];
%Return states as [z;S]
zdot = [ode(t,z(1:dopts.n),theta);S(:)];


function status = odeMaxTime(~,~,flag,maxTime)
%Callback to stop ODE if max time exceeded
persistent tt
status=0; %continue
switch(flag)
    case 'init'
        tt = tic;
    case ''
        if(toc(tt) > maxTime)
            fprintf(2,'ODE Integrator Max Time [%gs] Exceeded\n',maxTime);
            status = 1;
        end
end



% OLD CVODES CODE
% %Setup CVODES if selected
% if(strcmpi(dopts.integrator,'cvodes'))
%     if(dopts.nz0)
%         error('CVODES [SUNDIALS] does not currently solve problems with unknown initial conditions');        
%     end
%     %Setup CVODES Problem
%     data.theta = x0(1:dopts.np);
%     data.odeFun = prob.ode;
%     options = CVodeSetOptions('UserData',data,'RelTol',1e-6);
%     CVodeInit(@CVodeFun, 'BDF', 'Newton', prob.xdata(1), z0, options);
%     FSAOptions = CVodeSensSetOptions('method','Simultaneous',...
%                                      'RelTol',1e-8,...
%                                      'ErrControl', true,...
%                                      'ParamField', 'theta',...
%                                      'ParamScales', data.theta);
%     CVodeSensInit(dopts.np, [], zeros(dopts.n,dopts.np), FSAOptions);
% end
% 
%     case 'cvodes'       
%         %Reinitialize CVODE
%         data.theta = theta(1:dopts.np);
%         data.odeFun = odeFun;
%         opts = CVodeSetOptions('UserData',data,'RelTol',1e-6);
%         CVodeReInit(tm(1),z0,opts);
%         FSAopts = CVodeSensSetOptions('method','Simultaneous',...
%                                      'RelTol',1e-8,...
%                                      'ErrControl', true,...
%                                      'ParamField', 'theta',...
%                                      'ParamScales', data.theta);
%         CVodeSensReInit(S0,FSAopts);
%         %Call CVODEs
%         [~,~,y,yS] = CVode(tm(2:end),'Normal');
%         if(dopts.reqGrad)
%             f = [zeros(1,dopts.n*dopts.np);reshape(yS,[dopts.n*dopts.np length(tm)-1])'];
%         else
%             f = [z0';y'];
%         end
%         indGrad = false;
% 
% function [zd, flag, new_data] = CVodeFun(t, z, data)
% %CVODES Callback Function
% zd = data.odeFun(t,z,data.theta);
% flag = 0;
% new_data = [];

