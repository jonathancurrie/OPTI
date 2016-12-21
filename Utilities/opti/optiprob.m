function prob = optiprob(varargin)
%OPTIPROB  Create or alter a Problem for Optimization with OPTI
%
% problem = optiprob('param1',value1,'param2',value2,...) creates an OPTI
% problem structure with the parameters 'param' set to their corresponding
% values in 'value'. Parameters not specified will be set to the OPTI
% default.
%
% problem = optiprob(oldprob,'param1',value1,...) creates a copy of the old
% problem 'oldprob' and then fills in (or writes over) the parameters
% specified by 'param' and 'value'.
%
% problem = optiprob() creates a problem structure with all fields set to
% OPTI defaults.
%
% optiprob() prints a list of all possible fields and their function.

%   Copyright (C) 2011-2013 Jonathan Currie (IPL)

%If empty print all possible arguments
if (nargin == 0) && (nargout == 0)
    printfields();
    return
end

%Defaults
Names = {'Name';'f';'H';'Hstr';'fun';'sense';'objbias';'A';'b';'Aeq';'beq';'rl';'ru';'lb';'ub';'Q';'l';'qrl';'qru';'sdcone';...
         'nlcon';'nljac';'nljacstr';'nlrhs';'nle';'cl';'cu';'int';'sos';'xdata';'ydata';'weighting';'x0';'probtype';'solver';'path';'opts';...
         'ode';'odez0'};
Defaults = {'OPTI Problem';[];[];[];[];1;0.0;[];[];[];[];[];[];[];[];[];[];[];[];[];[];[];[];[];[];[];[];[];[];[];[];[];[];[];[];[];[];[];[]};         
altPairs = {{'grad','f'},{'c','f'},{'obj','fun'},{'hess','H'},{'ineq','A','b'},{'lin','Ar','rl','ru'},{'ndec','x0'},...
            {'jac','nljac'},{'jacstr','nljacstr'},{'hessstr','Hstr'},...
            {'eq','Aeq','beq'},{'bounds','lb','ub'},{'mix','A','b','e'},{'qp','H','f'},{'sdp','sdcone'},{'sedumi','sdcone'},...
            {'ctype','int'},{'vtype','int'},{'xtype','int'},{'ivars','int'},{'nlmix','nlcon','nlrhs','nle'},{'nl','nlcon','cl','cu'},...
            {'qc','q','l','qru'},{'qcrow','q','l','qrl','qru'},{'sle','A','b'},{'data','xdata','ydata'},...
            {'weights','weighting'},{'options','opts'},{'problemtype','probtype'},{'odefun','ode'},{'z0','odez0'},{'theta0','x0'},{'nleq','nlcon'}};

%Collect Sizes and lowercase matches         
m = size(Names,1); numberargs = nargin;
%Create structure with all names and default values
st = [Names,Defaults]'; prob = struct(st{:});

try
    % Check we have char or structure input. If structure, insert arguments
    i = 1;
    while i <= numberargs
        arg = varargin{i};
        if ischar(arg)
            break;
        end
        if ~isempty(arg)
            if ~isa(arg,'struct')
                error('An argument was supplied that wasn''t a string or struct!');
            end
            for j = 1:m
                if any(strcmp(fieldnames(arg),Names{j,:}))
                    val = arg.(Names{j,:});
                else
                    val = [];
                end
                if ~isempty(val)
                    checkfield(Names{j,:},val);
                    prob.(Names{j,:}) = val;
                end
            end
        end
        i = i + 1;
    end

    %Go through each argument pair and assign to correct field
    expectval = 0; %first arg is a name
    while i <= numberargs
        arg = varargin{i};
        switch(expectval)
            case 0 %field
                if ~ischar(arg)
                    error('Expected field name to be a string! - Argument %d',i);
                end
                oldfield = arg;
                %Check for special cases
                switch(lower(arg))
                    case 'ndec'
                        if(isempty(prob.x0))
                            prob.x0 = NaN(varargin{i+1},1); %don't need this if we have x0! [bugfix 20/6/13]
                        end
                        i = i + 1;
                        expectval = 0;
                    case 'sos'
                        expectval = 3; %maintain backwards compatibility by expecting 3 args (type, ind, wt)                    
                    otherwise
                        %Normal code to find match            
                        j = find(strcmpi(arg,Names) == 1);
                        if isempty(j)  % if no matches
                            %See if we have an alternative pair name
                            [j,mode] = checkPair(arg,Names,altPairs);
                            switch(mode) %switch on number in alt pair list
                                case 2
                                    if(strcmpi(arg,'sedumi') && ~isstruct(varargin{i+1})) %sedumi individual args
                                        expectval = 40;
                                    else
                                        expectval = 1; %just alternative name
                                    end
                                case 3 %e.g. qp,H,f
                                    if(strcmpi(arg,'qp')) %qp has H, f reversed (ahh hindsight..)
                                        expectval = -2;
                                    else
                                        expectval = 2; %two args eg ineq, eq
                                    end
                                case {4,5} %e.g. lin,A,rl,ru or qcrow,Q,l,qrl,qru
                                    expectval = mode-1;
                                otherwise
                                    error('Unknown OPTI problem form');
                            end
                        else
                            expectval = 1; % next arg is a value
                        end 
                end
            case 1 %val
                if(ischar(arg)), arg = lower(arg); end
                checkfield(Names{j,:},arg);            
                prob.(Names{j,:}) = arg;
                expectval = 0;
            case 2 %val 2x
                i = i + 1; arg2 = varargin{i};
                if(ischar(arg)), arg = lower(arg); end                
                if(ischar(arg2)), arg2 = lower(arg2); end                
                checkfield(Names{j,:},arg);
                prob.(Names{j,:}) = arg;                
                checkfield(Names{j+1,:},arg2);
                prob.(Names{j+1,:}) = arg2;
                expectval = 0;  
            case -2 %val 2x backwards
                i = i + 1; arg2 = varargin{i};
                if(ischar(arg)), arg = lower(arg); end
                if(ischar(arg2)), arg2 = lower(arg2); end                 
                checkfield(Names{j,:},lower(arg));
                prob.(Names{j,:}) = lower(arg);                
                checkfield(Names{j-1,:},arg2);
                prob.(Names{j-1,:}) = arg2;
                expectval = 0;
            case 3 %pre process, nonlinear mix, or sos
                switch(lower(oldfield))
                    case 'nlmix'            	
                        nlcon = arg; nlrhs = varargin{i+1}; nle = varargin{i+2};
                        checkfield('nlcon',nlcon); checkfield('nlrhs',nlrhs);
                        checkfield('nle',nle);
                        prob.nlcon = nlcon; prob.nlrhs = nlrhs; prob.nle = nle;
                        i = i + 2;
                        expectval = 0;
                    case 'nl'
                        nlcon = arg; cl = varargin{i+1}; cu = varargin{i+2};
                        checkfield('nlcon',nlcon); checkfield('cl',cl); checkfield('cu',cu);
                        prob.nlcon = nlcon; prob.cl = cl; prob.cu = cu;
                        i = i + 2;
                        expectval = 0;
                    case 'mix'
                        [A,b,Aeq,beq] = mix2gen(arg,varargin{i+1},varargin{i+2});
                        checkfield('A',A); checkfield('b',b);
                        checkfield('Aeq',Aeq); checkfield('beq',beq);
                        prob.A = A; prob.b = b;
                        prob.Aeq = Aeq; prob.beq = beq;
                        i = i + 2;
                        expectval = 0;
                    case 'lin'
                        A = arg; rl = varargin{i+1}; ru = varargin{i+2};
                        checkfield('A',A); checkfield('rl',rl); checkfield('ru',ru);
                        prob.A = A; prob.rl = rl; prob.ru = ru;
                        i = i + 2;
                        expectval = 0;
                     case 'sos'
                        if(isempty(arg))
                            %Check if we have 1 or 3 inputs
                            if(numberargs >= i + 2) %means 3 args OR not at end
                                if(~ischar(varargin{i+1})) %has to be 3 args
                                    i = i + 2;
                                end
                            end %has to be 1 arg                        
                        else
                            if(isstruct(arg))
                                checkfield('sos',arg);
                                sostype = arg.type; sosind = arg.index; soswt = arg.weight; 
                            elseif(numberargs >= i + 2)
                                sostype = arg; sosind = varargin{i+1}; soswt = varargin{i+2};
                                i = i + 2;
                            else
                                error('Unknown form of SOS constraints. Expected a structure, or 3 arguments (type, index, weight)');
                            end
                            checkfield('sos type',sostype); checkfield('sos index',sosind); checkfield('sos weight',soswt);
                            prob.sos = struct('type',sostype,'index',{sosind},'weight',{soswt});       
                        end
                        expectval = 0;    
                    case 'qc'
                        Q = arg; l = varargin{i+1}; qru = varargin{i+2};                    
                        %Auto fill in lower quad bound as -inf
                        if(iscell(qru))
                            qrl = cell(size(qru));
                            for j = 1:length(qrl)
                                qrl{j} = -Inf;
                            end
                        else
                            qrl = -Inf(size(qru));
                        end
                        checkfield('q',Q); checkfield('l',l); 
                        checkfield('qrl',qrl); checkfield('qru',qru);
                        prob.Q = Q; prob.l = l; prob.qrl = qrl; prob.qru = qru;     
                        i = i + 2;
                        expectval = 0;
                end
            case 4 %qc
                switch(lower(oldfield))
                    case 'qcrow'
                        %Expect all 4 args
                        Q = arg; l = varargin{i+1};
                        qrl = varargin{i+2};
                        qru = varargin{i+3};                                    
                        checkfield('q',Q); checkfield('l',l); 
                        checkfield('qrl',qrl); checkfield('qru',qru);
                        prob.Q = Q; prob.l = l; prob.qrl = qrl; prob.qru = qru;   
                        i = i+3;  
                        expectval = 0;
                end

            %SEDUMI Individual Args
            case 40
                %Get & Check Args
                At = arg; b = varargin{i+1}; c = varargin{i+2}; k = varargin{i+3};
                checkfield('SeDuMi At',At); checkfield('SeDuMi b',b); checkfield('SeDuMi c',c); checkfield('SeDuMi K',k);
                if(~isempty(prob.sdcone))
                    error('Cannot assign SeDuMi structure to field sdcone as it already contains data!');
                else
                    prob.sdcone = struct('At',At,'b',b,'c',c,'K',k);
                end
                i = i +3;
                expectval = 0;

        end
        i = i + 1;
    end

    if expectval %fallen off end somehow
        try
            error('Missing value for %s',arg);
        catch %#ok<CTCH>
            display(arg);
            error('Missing Field / Value for the above!');
        end
    end
catch ME
    throw(ME);
end


function [j,mode] = checkPair(arg,names,pairs)
%Check if user has supplied an alternative field name
no = length(pairs);
for i = 1:no
    if(strcmpi(arg,pairs{i}{1}))
        j = find(strcmp(pairs{i}{2},names) == 1);
        mode = length(pairs{i});
        return;
    end
end
throw(MException('OPTI:OPTIPROB_CHECKPAIR','Unrecognised parameter: %s',arg));


function checkfield(field,value)
%Check a field contains correct data type and values
if isempty(value), return; end %empty matrix is always valid

switch lower(field)
    %Function Handle or Matrix
    case {'f','h'}
        err = opticheckval.checkNumorFHndl(value,field);
    %Function Handle
    case {'fun','hmul','hstr','nlcon','nljac','nljacstr','ode'}
        err = opticheckval.checkFunHandle(value,field);
    %-1 or 1
    case 'sense'
        err = opticheckval.checkScalarSet(value,field,[-1 1]);
    %Double Scalar
    case 'objbias'
        err = opticheckval.checkScalarDbl(value,field);
    %Double Vector
    case {'b','beq','rl','ru','lb','ub','nlrhs','nle','cl','cu','x0','sedumi b','odez0'} 
        err = opticheckval.checkDblVec(value,field);
    %Double Vector (but may be sparse)
    case 'sedumi c'
        err = opticheckval.checkDblSpVec(value,field);
    %Double Matrix
    case {'a','aeq','sedumi at'} 
        err = opticheckval.checkDblMat(value,field);        
    %Double Vector or Cell array of Double Vectors
    case {'qrl','qru','sos index','sos weight'}
        err = opticheckval.checkDblVecOrCell(value,field);
    %Double Matrix or Cell array of Double Matrices
    case {'q','l','xdata','ydata','weighting'}
        err = opticheckval.checkDblMatOrCell(value,field);
    %SOS special struct
    case 'sos'
        err = opticheckval.checkStructFields(value,field,{'type','index','weight'});
    %SDCONE special struct
    case 'sdcone'
        if(isnumeric(value))
            err = opticheckval.checkDblMat(value,field);
        elseif(iscell(value))
            err = opticheckval.checkDblMatOrCell(value,field);
        elseif(isstruct(value))
            err = opticheckval.checkStructFields(value,field,{'At','b','c','K'});
        else
            err = MException('OPTI:SetFieldError','Parameter ''%s'' should be double matrix, cell array of double matrices, or SeDuMi input structure (At,b,c,K)',field);
        end    
    %Char or Numeric for Int
    case 'int'
        if(isnumeric(value))
            err = opticheckval.checkVectorIntGrtZ(value,field);
        elseif(ischar(value))
            value = lower(value);
            if(~all(value == 'c' | value == 'b' | value == 'i'))
            	err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a char array (string) with the following characters only:\n''c''\n''i''\n''b''\n',field);
            else
                err = [];
            end
        else
            err = MException('OPTI:SetFieldError','Parameter ''%s'' should be double vector of indices or char array of variable types',field);
        end                 
    %Char Array
    case {'name','sos type','path'}
        err = opticheckval.checkChar(value,field);
    %Structure
    case {'opts','sedumi k'}
        err = opticheckval.checkStruct(value,field);
    %Misc Strings
    case 'probtype'
        err = opticheckval.checkValidString(value,field,optiSolver('ptypes'));
    case 'solver'
        err = opticheckval.checkValidString(value,field,optiSolver('all'));
    otherwise  
        err = MException('OPTI:SetFieldError','Unrecognized parameter name ''%s''.', field);
end
if(~isempty(err)), throw(err); end


function printfields()
%Print out fields with answers
fprintf('\n OBJECTIVE FIELDS:\n');
fprintf('             fun: [ Nonlinear Objective Function (NL) ] \n');
fprintf('               f: [ LP/QP Linear Vector OR Gradient Function (NL) ] \n');
fprintf('               H: [ QP Quadratic Matrix OR Hessian of the Lagrangian Function (NL) ] \n');
fprintf('            Hstr: [ Hessian Function Structure (NL) ] \n');
fprintf('           sense: [ Minimization (1) or Maximization (-1) ] \n');
fprintf('         objbias: [ Linear or Quadratic Objective Bias Term ] \n');

fprintf('\n LINEAR CONSTRAINT FIELDS:\n');
fprintf('               A: [ Linear Inequality LHS (A*x <= b) ] \n');
fprintf('               b: [ Linear Inequality RHS ] \n');
fprintf('             Aeq: [ Linear Equality LHS (Aeq*x = beq) ] \n');
fprintf('             beq: [ Linear Equality RHS ] \n');
fprintf('              rl: [ Linear Constraint Lower Bound (rl <= A*x <= ru) ] \n');
fprintf('              ru: [ Linear Constraint Upper Bound ] \n');

fprintf('\n NONLINEAR CONSTRAINT FIELDS:\n');
fprintf('           nlcon: [ Nonlinear Constraint Function (NL) ] \n');
fprintf('           nlrhs: [ Nonlinear Constraint RHS (NL) ] \n');
fprintf('             nle: [ Nonlinear Constraint Types (-1, 0, 1) (NL) ] \n');
fprintf('              cl: [ Nonlinear Constraint Lower Bound (cl <= nlcon(x) <= cu) ] \n');
fprintf('              cu: [ Nonlinear Constraint Upper Bound ] \n');
fprintf('           nljac: [ Nonlinear Constraint Jacobian (NL) ] \n');
fprintf('        nljacstr: [ Nonlinear Constraint Jacobian Structure (NL) ] \n');
fprintf('               Q: [ Quadratic Constraint Quadratic LHS (x''Qx + l''x <= qru) ] \n');
fprintf('               l: [ Quadratic Constraint Linear LHS ] \n');
fprintf('             qrl: [ Quadratic Constraint Lower Bound (qrl <= x''Qx + l''x <= qru) ] \n');
fprintf('             qru: [ Quadratic Constraint Upper Bound ] \n');

fprintf('\n OTHER CONSTRAINT FIELDS:\n');
fprintf('              lb: [ Lower Bounds (lb <= x <= ub) ] \n');
fprintf('              ub: [ Upper Bounds ] \n');
fprintf('             int: [ Binary/Integer Variable String (''C'', ''I'' or ''B'')] \n');
fprintf('             sos: [ SOS Constraint Structure (.type .index .weight) ] \n');
fprintf('          sdcone: [ Semidefinite Cone Constraint (F1*x1 + F2*x2 + ... + Fn*xn - F0 >= 0 [PSD]) ] \n');

fprintf('\n DATA FITTING FIELDS:\n');
fprintf('           xdata: [ Data Fitting Problem X Data ] \n');
fprintf('           ydata: [ Data Fitting Problem Y Data ] \n');
fprintf('         weights: [ Data Fitting Problem Weights on Y Data ] \n');

fprintf('\n OTHER FIELDS:\n');
fprintf('            Name: [ Problem Name ] \n');
fprintf('              x0: [ Initial Solution Guess ] \n');
fprintf('        probtype: [ Optional string containing the problem type to be solved ] \n');

fprintf('\n DYNAMIC OPTIMIZATION FIELDS:\n');
fprintf('             ode: [ ODE Function ] \n');
fprintf('           odez0: [ ODE Initial Conditions ] \n');

%Print out groupings
fprintf('\n OBJECTIVE GROUPS:\n')
fprintf('          Quadratic Program:    ''qp''      [ opti(''qp'',H,f) ]\n');
fprintf(' System of Linear Equations:    ''sle''     [ opti(''sle'',A,b) ]\n');

fprintf('\n CONSTRAINT GROUPS:\n')
fprintf('   Decision Variable Bounds:    ''bounds''  [ opti(''bounds'',lb,ub) ]\n');
fprintf('        Linear Inequalities:    ''ineq''    [ opti(''ineq'',A,b) ]\n');
fprintf('          Linear Equalities:    ''eq''      [ opti(''eq'',Aeq,beq) ]\n');
fprintf('     Linear Row Constraints:    ''lin''     [ opti(''lin'',A,rl,ru) ]\n');
fprintf('               Linear Mixed:    ''mix''     [ opti(''mix'',A,b,e) ]\n');
fprintf('            Nonlinear Mixed:    ''nlmix''   [ opti(''nlmix'',nlcon,nlrhs,nle) ]\n');
fprintf('  Nonlinear Row Constraints:    ''nl''      [ opti(''nl'',nlcon,cl,cu) ]\n');
fprintf('     Quadratic Inequalities:    ''qc''      [ opti(''qc'',Q,l,qru) ]\n');
fprintf('  Quadratic Row Constraints:    ''qcrow''   [ opti(''qcrow'',Q,l,qrl,qru) ]\n');
fprintf('  Special Ordered Set (SOS):    ''sos''     [ opti(''sos'',type,index,weight) ]\n');
fprintf('      Data Fitting Data Set:    ''data''    [ opti(''data'',xdata,ydata) ]\n');
fprintf('           SeDuMi Arguments:    ''sedumi''  [ opti(''sedumi'',At,b,c,K) ]\n');
