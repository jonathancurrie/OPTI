function cb = buildCFun(mode,sobj,svar,opts,ncon)
% BUILDCFUN Build a C Code function for nonlinear callbacks

if(nargin > 3 && ~isempty(opts) && isfield(opts,'cbmode') && strcmpi(opts.cbmode,'cppad'))
    tname = 'symb_cadtemp.cpp';
    fname = 'symb_ccb.cpp';
    cmode = 'A';
else
    tname = 'symb_ctemp.c';
    fname = 'symb_ccb.c';
    cmode = 'C';
end
cb = [];

%Determine callback type
switch(mode)
    case 'new'
        %If compiling against CppAD, must use Intel or Visual Studio compiler
        if(cmode=='A')
            SymBuilder.CheckCppADCompile(true);
        elseif(cmode=='C') %check for LCC, doesn't work (don't know why)
            c = mex.getCompilerConfigurations('c');
            if(~isempty(c) && ~isempty(strfind(c.Name,'Lcc')))
                error('Lcc does not correctly compile SymBuilder C-Code callbacks. Please use Visual Studio or Windows SDK.');
            end
        end            
        %Find template
        p = which(tname);
        if(isempty(p)), error('Cannot find C template'); end
        %Copy template
        str = fileread(p);
        s2 = sprintf('Symbolic Builder Auto-Generated Callback Function\n');
        s2 = sprintf('%s * Generated %s',s2,datestr(now));
        str = regexprep(str,'SYMB_CTEMP - Template for generating SymBuilder C Code Callbacks',s2);
        fp = fopen(fname,'w');
        fprintf(fp,'%s\n\n',str);
        fclose(fp);
        return;
    case 'obj'
        title = 'Objective Function Callback';
        if(cmode=='A')
            fcall = sprintf('template <typename Type> Type objective(const vector<Type> &x)');
        else
            fcall = sprintf('double objective(double *x)');
        end
        var = 'j';
        cb = @(x) symb_ccb('obj',x);
        if(opts.verbose), fprintf('Generating Objective....'); end
    case 'grad'
        title = 'Objective Gradient Callback';
        fcall = sprintf('void gradient(double *x, double *v)');
        var = 'g';
        cb = @(x) symb_ccb('grad',x);
        if(nargin > 1 && ~isempty(sobj))
            if(opts.verbose), fprintf('Generating Gradient....'); end
        end
    case 'con'
        title = 'Constraint Function Callback';
        if(cmode=='A')
            fcall = sprintf('template <typename Type> void constraints(const vector<Type> &x, vector<Type> &v)');
        else
            fcall = sprintf('void constraints(double *x, double *v)');
        end            
        var = 'c';
        cb = @(x) symb_ccb('con',x);
        if(nargin > 1 && ~isempty(sobj))
            if(opts.verbose), fprintf('Generating Constraints....'); end 
        end
    case 'jac'
        title = 'Constraint Jacobian Callback';
        fcall = sprintf('void jacobian(double *x, double *pr, mwIndex *ir, mwIndex *jc)');
        var = 'J';
        cb = @(x) symb_ccb('jac',x);
        if(nargin > 1 && ~isempty(sobj))
            if(opts.verbose), fprintf('Generating Jacobian....'); end
        end
    case 'hess'
        title = 'Hessian of the Lagrangian Callback';
        fcall = sprintf('void hessian(double *x, double sigma, double *lambda, double *pr, mwIndex *ir, mwIndex *jc)');
        var = 'H';
        cb = @(x,sigma,lambda) symb_ccb('hess',x,sigma,lambda);
       if(nargin > 1 && ~isempty(sobj))
           if(opts.verbose), fprintf('Generating Hessian....'); end
       end
    otherwise
        error('Unknown MFun mode: %s',mode);
end

%Open file
fp = fopen(fname,'a+');
%Write Header   
fprintf(fp,'// %s\n',title);
fprintf(fp,fcall);
fprintf(fp,'\n{\n');

%If we have a function to write
if(nargin > 1 && ~isempty(sobj))
    %Need to rename all variables to x[0]...x[n]
    xvar = cell(length(svar),1);
    for i = 1:length(xvar)
        xvar{i} = sprintf('x[%d]',i);
    end     
    %Ensure both columns
    if(size(svar,1) > 1), svar = svar.'; end
    if(size(xvar,1) > 1), xvar = xvar'; end
    %Subs out individual symbolic variables into our indexed list and converts to normal numbers
    wstate = warning('off','symbolic:sym:sym:DeprecateExpressions');
    eq = vpa(subs(sobj,svar,xvar),16); %this takes too long - any suggestions?
    warning(wstate);
    %Sub out Lambda if Hessian
    if(var=='H')
        l = cell(ncon,1);
        l2 = cell(ncon,1);
        for i = 1:ncon
            l{i} = sprintf('lambda(%d)',i);
            l2{i} = sprintf('lambda[%d]',i);
        end
        wstate = warning('off','symbolic:sym:sym:DeprecateExpressions');
        eq = subs(eq,l,l2);
        warning(wstate);
    end

    %Enter Equations (Var Type Dictates Entry Type)
    switch(var)
        case 'j' %scalar
            fprintf(fp,'%s\n',ppc(cmode,regexprep(ccode(eq),'t0 = ','return ')));
            %Spaces at end
            fprintf(fp,'}\n\n');
            %Extra #Var function
            writeVarFcn(fp,length(svar));
        case {'g','c'} %vector, dense    
            %HACK
            if(~isa(eq,'sym')), eq = sym(eq); end
            if(size(eq,2) > 1), eq = eq.'; end
%             %Convert equations to string, write each as we go
%             for i = 1:length(eq)
%                 fprintf(fp,'%s\n',ppc(cmode,regexprep(ccode(eq(i)),'t0 = ',sprintf('v[%d] = ',i-1))));
%             end  
            if(length(eq) > 1) %write all equations at once
                fprintf(fp,'%s\n',ppc(cmode,regexprep(ccode(eq),'eq[(\d*)\][0\]','v[$1\]')));
            else
                fprintf(fp,'%s\n',ppc(cmode,regexprep(ccode(eq),'t0 = ','v[0] = ')));
            end
            if(var=='c')
                %Spaces at end
                fprintf(fp,'}\n\n');
                %Extra #Con function
                writeConFcn(fp,length(eq));
            end
        case {'J','H'} %matrix, sparse
            %Get NonZero elements from equation
            nzel = logical(eq ~= 0);               
            [rows,cols] = find(sparse(nzel));
            nz = nnz(nzel);                
            nzsym = eq(nzel);
            %Write Ir (rows)
            for i = 1:length(rows)
                fprintf(fp,'  ir[%d] = %d;\n',i-1,rows(i)-1);
            end
            %Write Jc (col starts)
            s = 0;
            for i = 1:size(eq,2)+1                
                fprintf(fp,'  jc[%d] = %d;\n',i-1,s);
                s = s + sum(cols==i);
            end            
            %Write Pr (vals/eqs)
%             for i = 1:length(nzsym)
%                 fprintf(fp,'%s\n',regexprep(ccode(nzsym(i)),'t0 = ',sprintf('pr[%d] = ',i-1)));
%             end
            if(length(nzsym) > 1)
                str = regexprep(ccode(nzsym),'nzsym[(\d*)\][0\]','pr[$1\]');
                if ~isempty(strfind(str, 'nzsym'))
                    nzsym = nzsym.';
                    str = regexprep(ccode(nzsym),'nzsym[(\d*)\][0\]','pr[$1\]');
                end
                fprintf(fp,'%s\n',ppc(cmode,str));
            else
                fprintf(fp,'%s\n',regexprep(ccode(nzsym),'t0 = ','pr[0] = '));
            end
            %Write extra NNZ function
            fprintf(fp,'}\n\n');
            writeNNZFcn(fp,nz,mode);
    end
else
    fprintf(fp,'   mexErrMsgTxt("%s is not supported in this SymBuilder Callback");\n',mode);
end
%Spaces at end
fprintf(fp,'}\n\n');
%NNZ default functions
if((nargin < 2 || isempty(sobj)) && (strcmpi(mode,'jac') || strcmpi(mode,'hess')))
    writeNNZFcn(fp,0,mode);
    fprintf(fp,'}\n\n');
end
%Default con function
if((nargin < 2 || isempty(sobj)) && strcmpi(mode,'con'))
    writeConFcn(fp,0);
    fprintf(fp,'}\n\n');
end
%Close file
fclose(fp);
if(nargin > 1 && ~isempty(sobj))
    if(opts.verbose), fprintf('Done\n'); end
end


%Write getNoVarFcn
function writeVarFcn(fp,len)
fprintf(fp,'//Return number of variables\n');
fprintf(fp,'mwIndex getNoVar()\n{\n');
fprintf(fp,'   return %d;\n',len);

%Write getNoConFcn
function writeConFcn(fp,len)
fprintf(fp,'//Return number of constraints\n');
fprintf(fp,'mwIndex getNoCon()\n{\n');
fprintf(fp,'   return %d;\n',len);

%Write getNNZFcn
function writeNNZFcn(fp,nz,mode)
if(strcmpi(mode,'jac'))
    fprintf(fp,'//Return number of nonzeros in the Jacobian\n');
    fprintf(fp,'mwIndex getNNZJac()\n{\n');
else
    fprintf(fp,'//Return number of nonzeros in the Hessian\n');
    fprintf(fp,'mwIndex getNNZHess()\n{\n');
end
fprintf(fp,'   return %d;\n',nz);

%Preprocess String to make it compatible with CppAD
function str = ppc(mode,str)

if(strcmpi(mode,'C')), return; end
if(isempty(strfind(str,'.0'))), return; end

str = strrep(str,'.0)',')');
str = strrep(str,'.0,',',');
str = strrep(str,'.0;',';');
str = strrep(str,'.0+','+');
str = strrep(str,'.0-','-');
str = strrep(str,'.0*','*');
str = strrep(str,'.0/','/');
str = strrep(str,'.0 ',' ');




