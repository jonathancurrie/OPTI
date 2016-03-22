function cb = buildMFun(mode,sobj,svar,opts)
% BUILDMFUN Build a MATLAB function for nonlinear callbacks

fname = 'symb_cb.m';
cb = [];

%Determine callback type
switch(mode)
    case 'new'
        title = 'SymBuilder Callback';
        if(~strcmpi(opts.use2ndDerivs,'yes'))
            fcall = sprintf('function v = symb_cb(mode,x)\n');
        else
            fcall = sprintf('function v = symb_cb(mode,x,sigma,lambda)\n');
        end
    case 'obj'
        title = 'Objective Function Callback';
        fcall = sprintf('function j = symb_%s(x)\n',mode);
        var = 'j';
        cb = @(x) symb_cb('obj',x);
        if(opts.verbose), fprintf('Generating Objective....'); end
    case 'grad'
        title = 'Objective Gradient Callback';
        fcall = sprintf('function g = symb_%s(x)\n',mode);
        var = 'g';
        cb = @(x) symb_cb('grad',x);
        if(opts.verbose), fprintf('Generating Gradient....'); end
    case 'con'
        title = 'Constraint Function Callback';
        fcall = sprintf('function c = symb_%s(x)\n',mode);
        var = 'c';
        cb = @(x) symb_cb('con',x);
        if(opts.verbose), fprintf('Generating Constraints....'); end
    case 'jac'
        title = 'Constraint Jacobian Callback';
        fcall = sprintf('function J = symb_%s(x)\n',mode);
        var = 'J';
        cb = @(x) symb_cb('jac',x);
        if(opts.verbose), fprintf('Generating Jacobian....'); end
    case 'hess'
        title = 'Hessian of the Lagrangian Callback';
        fcall = sprintf('function H = symb_%s(x,sigma,lambda)\n',mode);
        var = 'H';
        cb = @(x,sigma,lambda) symb_cb('hess',x,sigma,lambda);
        if(opts.verbose), fprintf('Generating Hessian....'); end
    otherwise
        error('Unknown MFun mode: %s',mode);
end

%Open file
if(strcmpi(mode,'new'))
    fp = fopen(fname,'w');
else
    fp = fopen(fname,'a+');
end
%Write Header            
fprintf(fp,fcall);
% fprintf(fp,'%% %s\n',upper(name));
fprintf(fp,'%% %s\n',title);
if(strcmpi(mode,'new'))
    fprintf(fp,'%% Symbolic Builder Auto-Generated Callback Function\n');
    fprintf(fp,'%% Generated %s\n\n',datestr(now));
    %Switch yard for function callbacks
    fprintf(fp,'%% Switch Yard\n');
    fprintf(fp,'switch(mode)\n');
    fprintf(fp,'    case ''obj'', v = symb_obj(x);\n');
    if(strcmpi(opts.use1stDerivs,'yes'))
        fprintf(fp,'    case ''grad'', v = symb_grad(x);\n');
    end
    if(opts.havCon)
        fprintf(fp,'    case ''con'', v = symb_con(x);\n');
        if(strcmpi(opts.use1stDerivs,'yes'))
            fprintf(fp,'    case ''jac'', v = symb_jac(x);\n');
        end
    end
    if(strcmpi(opts.use2ndDerivs,'yes'))
        fprintf(fp,'    case ''hess'', v = symb_hess(x,sigma,lambda);\n');
    end
    fprintf(fp,'    otherwise\n        error(''Unknown callback mode'')\n');
    fprintf(fp,'end\n\n\n');
	fclose(fp);
    return;
end
    

%Distribute Vars            
v = SymBuilder.detVarNames(svar);
if(~strcmp(v{1,1},'x') || size(v,1) > 1)
    fprintf(fp,'%% Slice Input Array\n');
    for i = 1:size(v,1)
        if(v{i,2}(end)-v{i,2}(1) == 0)
            fprintf(fp,'%s = x(%d);\n',v{i,1},v{i,2}(1));
        else
            fprintf(fp,'%s = x(%d:%d);\n',v{i,1},v{i,2}(1),v{i,2}(end));
        end
    end
end

%Check for preallocation
if(nargin < 4 || ~isstruct(opts))
    preallocate = true;
elseif(isfield(opts,'preallocate') && strcmpi(opts.preallocate,'yes'))
    preallocate = true;
else
    preallocate = false;   
end

%Build cell array of indexed variable strings
ivar = buildVarIndex(v,length(svar));      
%Ensure both columns
if(size(svar,1) > 1), svar = svar.'; end
if(size(ivar,1) > 1), ivar = ivar'; end
%Subs out individual symbolic variables into our indexed list and converts to normal numbers
wstate = warning('off','symbolic:sym:sym:DeprecateExpressions');
eq = vpa(subs(sobj,svar,ivar),16);
warning(wstate);

%Get equation size (matrices treated differently)
[r,c] = size(eq);

%Enter Equations (Var Type Dictates Entry Type)
switch(var)
    case 'j' %scalar
        fprintf(fp,'\n%% Equation:\n');
        fprintf(fp,'%s(1,1) = %s;\n',var,char(eq));
    case {'g','c'} %vector, dense    
        %HACK
        if(~isa(eq,'sym')), eq = sym(eq); end
        %Check for row
        if(c > 1), eq = eq.'; tr = 1; else tr = 0; end
        %Convert equations to string
        str = char(eq);
        %Remove 'matrix' if is found
        ind = strfind(str,'matrix');
        if(ind)
            str = str(ind+7:end-1);
        end
        if(preallocate)
            fprintf(fp,'\n%% Preallocate\n');
            if(tr)
                fprintf(fp,'%s = zeros(1,%d);\n',var,length(eq));
            else
                fprintf(fp,'%s = zeros(%d,1);\n',var,length(eq));
            end
        end
        if(strcmp(str(1),'[')), str = str(2:end); end %remove extra square brackets
        if(strcmp(str(end-1:end),']]')), str = str(1:end-2); end
        ss = regexp(str,'], ','split');
        fprintf(fp,'\n%% Equations:\n');
        if(tr)
            for i = 1:length(ss)
                if(~strcmp(ss{i}(2:end),'0'))
                    if(strcmp(ss{i}(1),'[')), ss{i} = ss{i}(2:end); end %remove extra square brackets
                    fprintf(fp,'%s(1,%d) = %s;\n',var,i,ss{i});
                end
            end
        else
            for i = 1:length(ss)
                if(~strcmp(ss{i}(2:end),'0'))
                    if(strcmp(ss{i}(1),'[')), ss{i} = ss{i}(2:end); end %remove extra square brackets
                    fprintf(fp,'%s(%d,1) = %s;\n',var,i,ss{i});
                end
            end   
        end
    case {'J','H'} %matrix, sparse
        %Get NonZero elements from equation
        nzel = logical(eq ~= 0);               
        [rows,cols] = find(sparse(nzel));
        nz = nnz(nzel);                
        nzsym = eq(nzel);

        %OLD WAY (Inefficient?)
        if(preallocate)
            fprintf(fp,'\n%% Preallocate\n');
            fprintf(fp,'%s = spalloc(%d,%d,%d);\n',var,r,c,nz);
        end
        fprintf(fp,'\n%% Sparse Matrix:\n');
        for i = 1:nz
            fprintf(fp,'%s(%d,%d) = %s;\n',var,rows(i),cols(i),char(nzsym(i)));
        end

        %FASTER BUT LESS CLEAR
%         %Row & Col Indices
%         fprintf(fp,'\n%% Row & Col Indices\n');
%         fprintf(fp,'rows = [');
%         for i = 1:length(rows)
%             if(i == length(rows))
%                 fprintf(fp,'%d];\n',rows(i));
%             else
%                 fprintf(fp,'%d, ',rows(i));
%             end
%         end
%         fprintf(fp,'cols = [');
%         for i = 1:length(cols)
%             if(i == length(cols))
%                 fprintf(fp,'%d];\n',cols(i));
%             else
%                 fprintf(fp,'%d, ',cols(i));
%             end
%         end
%         fprintf(fp,'\n%% Preallocate\n');
%         fprintf(fp,'vals = zeros(%d,1);\n',nz);
%         
%         fprintf(fp,'\n%% Nonzero Entries\n');
%         for i = 1:nz
%             fprintf(fp,'vals(%d,1) = %s;\n',i,char(nzsym(i)));
%         end
%         
%         fprintf(fp,'\n%% Create Sparse Return Matrix\n');
%         fprintf(fp,'%s = sparse(rows,cols,vals);\n',var);
end
%Spaces at end
fprintf(fp,'\n\n');
%Close file
fclose(fp);
if(opts.verbose), fprintf('Done\n'); end
rehash; %a shame, but seems to be required


%Convert variables from names to indexed strings [nominally m1 to m(1)]
function ivar = buildVarIndex(v,no)   

%Create cell to store strings in
ivar = cell(no,1); n = 1;

%Process v cell matrix, building strings as we go
for i = 1:size(v,1)
    name = v{i,1};
    indices = v{i,2};
    %If single variable, easy!
    if(length(indices) == 1)
        ivar{n} = name; n = n + 1;
    %Otherwise have to run through indexing list
    else
        for j = 1:length(indices)
            ivar{n} = sprintf('%s(%d)',name,j); n = n + 1;
        end
    end
end      