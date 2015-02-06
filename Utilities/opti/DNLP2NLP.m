function prob = DNLP2NLP(prob,opts)
%DNLS2NLS  Converts a Dynamic NLP problem into a NLP problem
% prob= DNLP2NLP(prob,opts)

%   Copyright (C) 2014 Jonathan Currie (I2C2)

% if(~isfield(prob,'type') || isempty(prob.type))
%     error('This function is not for user use, and should only be called from OPTI');
% end

%Get warning level
% warn = optiWarnLevel(opts.warnings);

%Get Derivative Checker Level
derCheck = false;
% if(strcmpi(opts.derivCheck,'on')), derCheck = true; end

%Transpose as required
if(size(prob.odez0,2) > 1), prob.odez0 = prob.odez0'; end

%Get Options & Set Sizes
dopts = optidynset(opts.dynamicOpts);
dopts.n = length(prob.odez0);           %number of states
dopts.nu = length(prob.x0);             %number of input steps (+ parameters?)

% %Set Maximum Integrator Time if no output function specified
% if(~isfield(dopts.odeOpts,'OutputFcn') || isempty(dopts.odeOpts.OutputFcn))
%     dopts.odeOpts = odeset(dopts.odeOpts,'OutputFcn',@(t,y,flag) odeMaxTime(t,y,flag,dopts.odeMaxTime));
% else
%     if(warn > 1), optiwarn('OPTI:NoOutFcn','ODE Output Function has already been specified by the user, ODE max time cannot be set'); end
% end
%Get Initial Time 
t0 = dopts.initialT;

%Check ODE Function given above sizes
try
    prob.ode(1,ones(dopts.n,1),ones(dopts.np,1));
catch ME
    throwAsCaller(MException('OPTI:ODE_ERROR','There was an error evaluating your ODE function. This is normally due to an incorrectly specified x0 and/or z0.\n\nCheck the following error for more information:\n\n%s',ME.message));
end

%Decide what sensitivity strategy to use
if(isempty(dopts.sensitivity))
    if(~isempty(dopts.dfdz) || ~isempty(dopts.dfdu))
        dopts.sensitivity = 'User';
    else
        dopts.sensitivity = 'ND';
    end
end

%Check User Derivatives (if supplied) and/or Setup Derivative Estimation Method if required
z0 = prob.odez0; u0 = prob.x0;
if(strcmpi(dopts.sensitivity,'user'))
    if(~isempty(dopts.dfdz))
        if(nargin(dopts.dfdz) ~= 3), error('ODE dfdz must be a function which accepts 3 arguments (t,z,u)'); end  
%         if(derCheck)
%             z = @(z) prob.ode(1,z,u0);
%             dfdz = @(z) dopts.dfdz(1,z,u0);
%             optiDerCheck(z,dfdz,z0,'ODE DFDZ',warn);
%         end
        dopts.dfdz_method = 2;
    else
%         if(warn>1), optiwarn('OPTI:NoUserDer','You have specified to use user supplied derivatives for DFDZ, but the corresponding function is empty.\nUsing ND instead.'); end
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
    if(~isempty(dopts.dfdu))
        if(nargin(dopts.dfdu) ~= 3), error('ODE dfdu must be a function which accepts 3 arguments (t,z,u)'); end
%         if(derCheck)
%             u = @(u) prob.ode(1,z0,u);
%             dfdu = @(u) dopts.dfdp(1,z0,u);
%             optiDerCheck(u,dfdu,u0,'ODE DFDU',warn);
%         end
        dopts.dfdu_method = 2;
    else
        if(warn>1), optiwarn('OPTI:NoUserDer','You have specified to use user supplied derivatives for DFDU, but the corresponding function is empty.\nUsing ND instead.'); end
        dopts.dfdu_method = 0;
    end
else
    %Setup derivative estimation method
    if(strcmpi(dopts.sensitivity,'cs'))
        dopts.dfdu_method = 3;
    elseif(strcmpi(dopts.sensitivity,'ad'))
        dopts.dfdu_method = 1;
    else
        dopts.dfdu_method = 0;
    end
end


%Assign Objective
dopts.reqGrad = false;
prob.fun = @(u,tm) odeEstim(prob.ode,prob.odez0,tm,u,dopts);
prob.misc.fitFun = prob.fun;

%Assign Gradient (if required)
if(~strcmpi(dopts.sensitivity,'none'))
    %Setup Initial Sensitivity
    dopts.Su0 = zeros(dopts.n,dopts.nu);
    
    %Add Derivative Calculation
    dopts.reqGrad = true;
    prob.f = @(u) odeEstim(prob.ode,prob.odez0,prob.xdata,u,dopts);
    prob.misc.fitGrad = prob.f;
end

function f = odeEstim(odeFun,z0,tm,u,dopts)
%Optimizer Callback Function   

%Check if we need gradient
if(dopts.reqGrad)
    ode = @(t,x) odeSens(t,x,odeFun,u(1:dopts.nu),dopts);        
    %Augment Sensitivity initial states
    Z0 = [z0;dopts.S0(:)];
else
    ode = @(t,x) odeFun(t,x,u(1:dopts.nu)); %original function
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
    if(r~=dopts.nT || c~=dopts.n*(1+dopts.nu))
        fprintf(2,'  INTEGRATOR ERROR IN GRADIENT EVALUATION - OPTI IS RETURNING 1e6 AT THIS POINT\n');
        f = 1e6*ones(dopts.nT,dopts.n*(1+dopts.nu));
    end
    %Index Sensitivity States (n+1) & States we are comparing to measured data
    f = f(:,dopts.inJacz);
    %Index Measurements within the Jacobian we will compare to
    if(dopts.selectedMeas)
        f = f(dopts.inJacMeas);
    end
    %Reshape into standard Jacobian
    f = reshape(f,dopts.nu,dopts.nu);    
else   
    %Check we got the correct sized vector, otherwise integrator may have failed
    [r,c] = size(f);
    if(r~=dopts.nT || c~=dopts.n)
        fprintf(2,'  INTEGRATOR ERROR IN OBJECTIVE EVALUATION - OPTI IS RETURNING 1e6 AT THIS POINT\n');
        f = 1e6*ones(dopts.nT,dopts.n);
    end
    %Run Objective Function
    f = dopts.obj(f);
end



function zdot = odeSens(t,z,ode,u,dopts)
%ODE with Integrated Sensitivity Function 
%Get dfdz
switch(dopts.dfdz_method)
    case 0 %numerical differentiation
        dfdz = mklJac(@(z) ode(t,z,u),z(1:dopts.n));
    case 1 %automatic differentiation
        dfdz = autoJac(@(z) ode(t,z,u),z(1:dopts.n));
    case 2 %user supplied
        dfdz = dopts.dfdz(t,z,u);
    case 3 %complex step
        dfdz = cstepJac(@(z) ode(t,z,u),z(1:dopts.n));
end
%Get dfdu
switch(dopts.dfdu_method)
    case 0 %numerical differentiation
        dfdu = mklJac(@(theta) ode(t,z(1:dopts.n),theta),u);
    case 1 %automatic differentiation
        dfdu = autoJac(@(theta) ode(t,z(1:dopts.n),theta),u);
    case 2 %user supplied
        dfdu = dopts.dfdp(t,z,u);
    case 3 %complex step
        dfdu = cstepJac(@(theta) ode(t,z(1:dopts.n),theta),u);
end
%Complete sensitivity differential equation
S = dfdz*reshape(z(dopts.n+1:end),dopts.n,dopts.nu) + dfdu;
%Return states as [z;S]
zdot = [ode(t,z(1:dopts.n),u);S(:)];


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

