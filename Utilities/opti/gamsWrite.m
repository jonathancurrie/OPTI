function gamsWrite(prob,filename)
%gamsWrite  Write a GAMS Model from MATLAB data
%
% gamsWrite(prob,filename) writes the OPTI problem to a GAMS model (.gms)
% and optionally a GDX file if SCIP is not used.
%
% You may specify a full path to the file, or if you specify a filename
% only, it will be written to the current directory.

%   Copyright (C) 2013 Jonathan Currie (Control Engineering)

ext = 'gms';

%Quick Checks
if(~ischar(filename))
    error('Filename must be a char array!');
end
if(~isempty(prob.sdcone))
    error('Semidefinite problems cannot be written by this interface');
end
% Make sure we have file extension for GAMS model
if(isempty(strfind(lower(filename),['.' ext])))
    gdxname = filename;
    filename = [filename '.' ext];    
else
    % Remove file extension for GDX file
    ind = strfind(lower(filename),['.' ext]);
    gdxname = filename(1:ind-1);
end 

% Check for full filepath
if(isempty(strfind(filename,':')))
    filename = [cd filesep filename];
    gdxname = [cd filesep gdxname];
end

%If we have SCIP installed, use that (much easier, and does nonlinear models)
if(optiSolver('scip',0))
    opts = optiset('solver','scip','warnings','none','display','off','solverOpts',scipset('gamsfile',filename));
    if(isstruct(prob.int) && isfield(prob.int,'str'))
        prob.int = prob.int.str;
    end
    if(isempty(prob.x0) || all(isnan(prob.x0))), prob.x0 = randn(prob.sizes.ndec,1); end
    Opt = opti(prob,opts);
    try
        solve(Opt);
    catch ME
        if(strcmpi(ME.identifier,'OPTI:SCIPVAR'))
            throwAsCaller(MException('OPTI:SCIPVAR','AN ERROR OCCURED WHILE PRE-PROCESSING YOUR MATLAB FUNCTION(S) - PLEASE EXAMINE THE SCIP ERROR BELOW:\n\n%s',ME.message));
        else
            throw(ME);
        end
    end
    return;
end

%Check not nonlinear
if(~isempty(prob.fun) && isa(prob.fun,'function_handle'))
    error('The OPTI GAMS Writer only works for problems up to MIQCQP complexity.\n\nFor general nonlinear problems you must have %s installed to write nonlinear GAMS problems','SCIP');
end

%Check GAMS routines are available
if(~exist(['rgdx.' mexext],'file') || ~exist(['wgdx.' mexext],'file'))
    error('Could not find the required GAMS routines on your PC [wgdx,rgdx]!\n\n%s',...
          'Please download a current GAMS distribution from "http://www.gams.com/download/", then add the GAMS install folder to the MATLAB path.');
end

%Check for SOS constraints (not yet implemented)
if(~isempty(prob.sos))
    optiwarn('OPTI:GamsWriteSOS','SOS constraints are not currently supported by this interface, so they will be skipped!');
end

%Determine number of decision variables
if(isfield(prob,'sizes'))
    ndec = prob.sizes.ndec;
elseif(~isempty(prob.f))
    ndec = numel(prob.f);
elseif(~isempty(prob.x0))
    ndec = numel(prob.x0);
elseif(~isempty(prob.H))
    ndec = size(prob.H,1);
else
    error('Could not determine the number of decision variables!');
end

%Convert linear constraints to row format
if(~isempty(prob.b) || ~isempty(prob.beq))
    [prob.A,prob.rl,prob.ru] = gen2row(prob.A,prob.b,prob.Aeq,prob.beq);
end

%Convert integer constraints to string
if(~isempty(prob.int) && ~isstruct(prob.int))
    if(ischar(prob.int))
        str = prob.int;
        prob.int = struct('str',str);
    elseif(isnumeric(prob.int))
        str = repmat('C',1,ndec);
        str(prob.int) = 'I';
        prob.int = struct('str',str);
    else
        error('Unknown integer format');
    end
elseif(isempty(prob.int))
    prob.int.str = [];
end
prob.int.str = upper(prob.int.str);
nint = sum(prob.int.str == 'I');
nbin = sum(prob.int.str == 'B');

%Starting writing GAMS Model File
fid = fopen(filename,'w');
if(fid < 0), error('Error creating GAMS Model: %s',filename); end

name = regexprep(prob.Name,' ','_');

fprintf(fid,'* OPTI Generated GAMS Model (%s)\n',datestr(now));
fprintf(fid,'*   Model Name             : %s\n',name);
if(isfield(prob,'type') && ~isempty(prob.type))
    fprintf(fid,'*   Model Classification   : %s\n',prob.type);
end
fprintf(fid,'*   Variables              : %d (%d Integer, %d Binary)\n',ndec,nint,nbin);
if(~isempty(prob.rl)),  fprintf(fid,'*   Linear Constraints     : %d\n',numel(prob.rl)); end
if(~isempty(prob.qrl)), fprintf(fid,'*   Quadratic Constraints  : %d\n',numel(prob.qrl)); end
if(~isempty(prob.cl)),  fprintf(fid,'*   Nonlinear Constraints  : %d\n',numel(prob.cl)); end
fprintf(fid,'\n');

try      
    model = 'lp';
    %Check for integer and binary constraints
    if(~all(prob.int.str == 'C'))
        int = struct('name','int','type','set','dim',1,'uels',{{find(prob.int.str ~= 'C')}},'form','full','ts','integer and binary indices');
        if(any(prob.int.str == 'B'))
            %Force binary bounds
            ind = (prob.int.str == 'B')';
            if(isempty(prob.lb)), prob.lb = -Inf(ndec,1); end
            if(isempty(prob.ub)), prob.ub = Inf(ndec,1); end
            %Set bin indicies which are infinite or bounds > 1 or < 0 to 1/0
            prob.lb((isinf(prob.lb) | (prob.lb < 0)) & ind) = 0;
            prob.ub((isinf(prob.ub) | (prob.ub > 1)) & ind) = 1;            
        end
    else
        int = [];
    end        
    %Check for bounds
    if(~isempty(prob.lb))
        lb = struct('name','lb','type','parameter','dim',1,'val',full(prob.lb),'form','full','ts','decision variable lower bound');        
    else
        lb = [];
    end
    if(~isempty(prob.ub))
        ub = struct('name','ub','type','parameter','dim',1,'val',full(prob.ub),'form','full','ts','decision variable upper bound');        
    else
        ub = [];
    end
    %Check for initial guess
    if(~isempty(prob.x0) && all(~isnan(prob.x0)))
        x0 = struct('name','x0','type','parameter','dim',1,'val',full(prob.x0),'form','full','ts','decision variable initial guess');
    else
        x0 = [];
    end
    %Check for linear objective
    if(~isempty(prob.f) && isnumeric(prob.f))
        f = struct('name','f','type','parameter','dim',1,'val',full(prob.f),'form','full','ts','linear objective vector');        
    else
        f = [];
    end
    %Check for quadratic objective
    if(~isempty(prob.H) && isnumeric(prob.H))
        %We require all elements of H, see if we have a saved copy before
        %triu/tril
        if(isfield(prob,'save') && isfield(prob.save,'H') && ~isempty(prob.save.H) && isnumeric(prob.save.H))
            H = prob.save.H;
        else
            H = prob.H;
        end
        [i,j,v] = find(sparse(H));
        H = struct('name','H','type','parameter','dim',2,'val',[i j v],'form','sparse','ts','quadratic objective matrix');
        model = 'qp';
    else
        H = [];
    end
    %Check for linear constraints
    if(~isempty(prob.A))   
        [i,j,v] = find(sparse(prob.A));
        A = struct('name','A','type','parameter','dim',2,'val',[i j v],'form','sparse','ts','linear constraint matrix');
        rl = struct('name','rl','type','parameter','dim',1,'val',full(prob.rl),'form','full','ts','linear constraint lower bound');
        ru = struct('name','ru','type','parameter','dim',1,'val',full(prob.ru),'form','full','ts','linear constraint upper bound');
    else
        A = []; rl = []; ru = [];
    end
    %Check for quadratic constraints
    if(~isempty(prob.Q))   
        %Concatenate each Q matrix into a 3D array
        if(iscell(prob.Q))
            i = []; j = []; k = []; v = [];
            for n = 1:length(prob.Q) %could preallocate...
                [I,J,V] = find(sparse(prob.Q{n}));
                i = [i;I]; %#ok<AGROW>
                j = [j;J]; %#ok<AGROW>
                k = [k;n*ones(size(I))]; %#ok<AGROW>
                v = [v;V]; %#ok<AGROW>
            end
        else
            [i,j,v] = find(sparse(prob.Q));
            k = ones(size(i));
        end
        Q = struct('name','Q','type','parameter','dim',3,'val',[i j k v],'form','sparse','ts','quadratic constraint matrix');
        [i,j,v] = find(sparse(prob.l));
        l = struct('name','l','type','parameter','dim',2,'val',[i j v],'form','sparse','ts','quadratic constraint linear matrix');
        qrl = struct('name','qrl','type','parameter','dim',1,'val',full(prob.qrl),'form','full','ts','quadratic constraint lower bound');
        qru = struct('name','qru','type','parameter','dim',1,'val',full(prob.qru),'form','full','ts','quadratic constraint upper bound');
        model = 'qcqp';
    else
        Q = []; l = []; qrl = []; qru = [];
    end

    %Write GAMS Data
    writeGDX(gdxname,H,f,lb,ub,x0,int,A,rl,ru,Q,l,qrl,qru);    

    %Write set definitions
    fprintf(fid,'Sets\n');
    fprintf(fid,' n /1*%d/\n',ndec);
    if(~isempty(int)), fprintf(fid,' ni(n)\n'); end
    if(~isempty(prob.rl)), fprintf(fid,' mlin /1*%d/\n',numel(prob.rl)); end
    if(~isempty(prob.qrl)), fprintf(fid,' mquad /1*%d/\n',numel(prob.qrl)); end
    if(~isempty(prob.cl)), fprintf(fid,' mnlin /1*%d/\n',numel(prob.cl)); end
    fprintf(fid,';\n');

    %Alias for QP, QC
    if(~isempty(H)), fprintf(fid,'Alias (n,n1);\n\n'); end

    %Write parameter definitions
    fprintf(fid,'Parameters\n');
    if(~isempty(H)), fprintf(fid,' H(n,n1)\n'); end
    if(~isempty(f)), fprintf(fid,' f(n)\n'); end
    if(~isempty(lb) && ~isempty(ub))
        fprintf(fid,' lb(n), ub(n)\n');
    elseif(~isempty(lb))
        fprintf(fid,' lb(n)\n');
    elseif(~isempty(ub))
        fprintf(fid,' ub(n)\n');
    end
    if(~isempty(x0)), fprintf(fid,' x0(n)\n'); end
    if(~isempty(A)), fprintf(fid,' A(mlin,n), rl(mlin), ru(mlin)\n'); end 
    if(~isempty(Q)), fprintf(fid,' Q(n,n1,mquad), l(n,mquad), qrl(mquad), qru(mquad)\n'); end
    fprintf(fid,';\n');

    %Enter Variables to Import
    fprintf(fid,'* Load Parameters from GDX\n');
    fprintf(fid,'$gdxin %s.gdx\n',gdxname);
    if(~isempty(H)), fprintf(fid,'$load H\n'); end
    if(~isempty(f)), fprintf(fid,'$load f\n'); end
    if(~isempty(lb) && ~isempty(ub))
        fprintf(fid,'$load lb,ub\n');
    elseif(~isempty(lb))
        fprintf(fid,'$load lb\n');
    elseif(~isempty(ub))
        fprintf(fid,'$load ub\n');
    end
    if(~isempty(int)), fprintf(fid,'$load ni = int\n'); end
    if(~isempty(x0)), fprintf(fid,'$load x0\n'); end
    if(~isempty(A)), fprintf(fid,'$load A,rl,ru\n'); end
    if(~isempty(Q)), fprintf(fid,'$load Q,l,qrl,qru\n'); end
    fprintf(fid,'$gdxin\n\n');

    %Complete GAMS file
    %Variable Definition
    fprintf(fid, 'Variables\n x(n)\n');
    fprintf(fid, ' obj;\n\n');

    %Integer and binary def
    if(~isempty(int))
        fprintf(fid,'Integer Variables\n x;\n');
        fprintf(fid,' x.prior(n) = inf;\n');
        fprintf(fid,' x.prior(ni) = 1;\n\n');
    end

    %Equations
    fprintf(fid, 'Equations\n');    
    if(~isempty(A))
        fprintf(fid, ' linconLB(mlin)\n'); 
        fprintf(fid, ' linconUB(mlin)\n'); 
    end
    if(~isempty(Q))
        fprintf(fid, ' quadconLB(mquad)\n');
        fprintf(fid, ' quadconUB(mquad)\n');
    end
    fprintf(fid, ' cost;\n\n');

    %Objective
    fprintf(fid, '* Cost Function\n');
    switch(model)
        case 'lp'
            fprintf(fid, ' cost.. obj =e= sum(n, f(n)*x(n));\n\n');
        case {'qp','qcqp'}
            fprintf(fid, ' cost.. obj =e= 0.5*sum(n,x(n)*sum(n1,H(n,n1)*x(n1))) + sum(n,f(n)*x(n));\n\n'); 
        otherwise
            error('Problem type ''%s'' not yet supported',prob.type);
    end

    %Bounds
    if(~isempty(lb) || ~isempty(ub))
        fprintf(fid, '* Decision Variable Bounds\n');
        if(~isempty(lb)), fprintf(fid, ' x.lo(n) = lb(n);\n'); end
        if(~isempty(ub)), fprintf(fid, ' x.up(n) = ub(n);\n'); end
        fprintf(fid,'\n');
    end

    %Initial Guess
    if(~isempty(x0))
        fprintf(fid, '* Initial Solution Guess\n');
        fprintf(fid, ' x.l(n) = x0(n);\n\n');
    end

    %Linear Constraints
    if(~isempty(A))
        fprintf(fid,'* Linear Constraints\n');     
        fprintf(fid,' linconLB(mlin).. rl(mlin) =l= sum(n, A(mlin,n)*x(n));\n');
        fprintf(fid,' linconUB(mlin).. sum(n, A(mlin,n)*x(n)) =l= ru(mlin);\n\n');            
    end

    %Quadratic Constraints
    if(~isempty(Q))
        fprintf(fid,'* Quadratic Constraints\n');
        fprintf(fid,' quadconLB(mquad).. qrl(mquad) =l= sum(n,x(n)*sum(n1,Q(n,n1,mquad)*x(n1))) + sum(n,l(n,mquad)*x(n));\n');
        fprintf(fid,' quadconUB(mquad).. sum(n,x(n)*sum(n1,Q(n,n1,mquad)*x(n1))) + sum(n,l(n,mquad)*x(n)) =l= qru(mquad);\n\n');        
    end
    
    %Write model declaration
    fprintf(fid,'\nModel %s / all /;\n\n',name);
    
    %Determine which GAMS solver type to use
    switch(model)
        case 'lp'
            if(isempty(int))
                fprintf(fid,'Solve %s using %s minimizing obj;\n',name,'LP');
            else        
                fprintf(fid,'Solve %s using %s minimizing obj;\n',name,'MIP');
            end
        case {'qp','qcqp'}
            if(isempty(int))
                fprintf(fid,'Solve %s using %s minimizing obj;\n',name,'QCP');
            else        
                fprintf(fid,'Solve %s using %s minimizing obj;\n',name,'MIQCP');
            end
    end

catch ME
    fclose(fid);
    rethrow(ME);
end

fclose(fid);
    

function writeGDX(name,H,f,lb,ub,x0,int,A,rl,ru,Q,l,qrl,qru)

args = {};
if(~isempty(H)),    args = [args H]; end
if(~isempty(f)),    args = [args f]; end
if(~isempty(lb)),   args = [args lb]; end
if(~isempty(ub)),   args = [args ub]; end
if(~isempty(x0)),   args = [args x0]; end
if(~isempty(int)),  args = [args int]; end
if(~isempty(A)),    args = [args A]; end
if(~isempty(rl)),   args = [args rl]; end
if(~isempty(ru)),   args = [args ru]; end
if(~isempty(Q)),    args = [args Q]; end
if(~isempty(l)),    args = [args l]; end
if(~isempty(qrl)),  args = [args qrl]; end
if(~isempty(qru)),  args = [args qru]; end

if(~isempty(args))
    wgdx(name,args{:});
end

