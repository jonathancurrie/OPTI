classdef scipvar
%SCIPVAR  SCIP Variable class
%
%   For internal SCIP MEX Interface use only. Allows conversion of a
%   deterministic MATLAB function into an SCIP MEX Interface compatible
%   instruction list. Use 'methods(scipvar)' to see available functions.
%
%   Copyright (C) 2012/2013 Jonathan Currie (Control Engineering)
    
    %#ok<*STOUT,*MANU,*INUSD>
    
    properties
        ins         % Instruction List
        indx        % Variable index
    end
    
    %Instructions
    properties (Constant)
       NUM = 0;
       VAR = 1;
       MUL = 3;
       DIV = 4;
       ADD = 5;
       SUB = 6;
       SQU = 7;
       SQR = 8;
       POW = 9;
       EXP = 10;
       LOG = 11;
       MIN = 15;
       MAX = 16;
       ABS = 17;
    end
    
    methods
        %SCIPVAR
        function s = scipvar(rows,cols)
        %SCIPVAR constructor
            switch(nargin)
                
                case 1
                    if(length(rows) == 2)
                        s = scipvar(rows(1),rows(2));
                    else
                        s = scipvar(rows,rows);
                    end
                    
                case 2
                    if(isnumeric(rows) && isnumeric(cols))
                        s(rows,cols) = scipvar; %preallocate
                        %Fill in indicies
                        k = 0;
                        for i = 1:cols
                            for j = 1:rows                            
                                s(j,i).indx = k; k = k + 1;
                            end
                        end
                    else
                        if(iscell(rows) && iscell(cols))
                            [r,c] = size(rows);
                            s(r,c) = scipvar; %preallocate
                            %Fill in indicies and existing instructions
                            for i = 1:r
                                for j = 1:c                 
                                    s(i,j).ins = rows{i,j};
                                    s(i,j).indx = cols{i,j};
                                end
                            end
                        else
                            error('Unknown constructor format');
                        end
                    end
            end
        end
        
        %-- SIMPLE OPERATORS --%
        % PLUS 
        function c = plus(a,b)
        %Implements vectorized a + b 
            c = vectorOp(a,b,scipvar.ADD);            
        end
        % MINUS 
        function c = minus(a,b)
        %Implements vectorized a - b 
            c = vectorOp(a,b,scipvar.SUB);            
        end
        % RDIVIDE 
        function c = rdivide(a,b)
        %Implements vectorized a ./ b 
            c = vectorOp(a,b,scipvar.DIV);            
        end
        % MRDIVIDE 
        function c = mrdivide(a,b)
        %Implements vectorized a ./ b if b is a scalar, otherwise error
            if(isscalar(b))
                c = vectorOp(a,b,scipvar.DIV);
            else
                error('This interface does not provide matrix right divide (A/b)');
            end
        end      
        % TIMES 
        function c = times(a,b)
        %Implements vectorized a .* b 
            c = vectorOp(a,b,scipvar.MUL);            
        end
        % MTIMES
        function c = mtimes(a,b)
        %Implements a*b or vectorized a.*b if a or b is scalar
            if(isscalar(a) || isscalar(b))
                c = vectorOp(a,b,scipvar.MUL);
            else
                c = matrixOp(a,b,scipvar.MUL);                    
            end
        end
        % POWER
        function c = power(a,b)
        %Implements a.^b
            c = vectorOp(a,b,scipvar.POW);
        end
        % MPOWER
        function c = mpower(a,b)
        %Implements A^b if b is an integer, a.^b if a is scalar, otherwise error
            if(isscalar(a))
                c = vectorOp(a,b,scipvar.POW);
            else
                if(isnumeric(b) && (round(b) == b)) %check for simple A*A*A...
                	[r,c] = size(a);
                    if(r == c)
                        c = a;
                        for i = 2:b
                            c = c * a;
                        end
                    else
                        error('You may only raise a square matrix to an integer power using this interface. Use .^ for element wise power.');
                    end                    
                else  
                    error('This interface only provides (A^b) when b is a integer scalar. Use .^ for element wise power.');
                end
            end
        end
        %UMINUS
        function c = uminus(a)
        %Implements -a
            c = -1*a;
        end
        
        %-- SIMPLE FUNCTIONS --%        
        % EXP
        function c = exp(a)
        %Implements vectorized exp(a) 
            c = vectorFcn(a,scipvar.EXP);
        end        
        % LOG
        function c = log(a)
        %Implements vectorized log(a)
            c = vectorFcn(a,scipvar.LOG);
        end        
        % LOG10
        function c = log10(a)
        %Implements vectorized log10(a)
            c = 0.4342944819032518.*vectorFcn(a,scipvar.LOG); %N.S. Identity
        end
        % ABS
        function c = abs(a)
        %Implements vectorized abs(a)
            c = vectorFcn(a,scipvar.ABS);
        end
        % SQRT
        function c = sqrt(a)
        %Implements vectorized sqrt(a)
            c = vectorFcn(a,scipvar.SQR);
        end
        % DOT 
        function c = dot(a,b)
        %Implements vectorized dot product of a, b
            c = sum(a.*b);
        end
        
        %-- EXPERIMENTAL FUNCTIONS --%        
        % NORM 
        function c = norm(a)
        %Implements vectorized norm(a) [2 norm for vectors, Frobenius norm for matrices)
            [row,col] = size(a);
            if(row > 1 && col > 1)
                c = sqrt(sum(sum(a.^2)));
            else
                c = sqrt(sum(a.^2));
            end
        end 
        
        %-- ARRAY FUNCTIONS --%
        % SUM
        function c = sum(a,n)
        %Implements sum(a)
            [r,c] = size(a);
            if(r > 1 || c > 1)
                if(c == 1) %sum down the rows
                    c = dimensionOp(a,1,scipvar.ADD);
                elseif(r == 1)%sum along the cols
                    c = dimensionOp(a,2,scipvar.ADD);
                else %use n as dimension
                    if(nargin == 1)
                        n = 1;
                    end
                    c = dimensionOp(a,n,scipvar.ADD);
                end
            else
                c = a;
            end
        end
        
        % PROD
        function c = prod(a,n)
        %Implements prod(a)
            [r,c] = size(a);
            if(r > 1 || c > 1)
                if(c == 1) %prod down the rows
                    c = dimensionOp(a,1,scipvar.MUL);
                elseif(r == 1)%prod along the cols
                    c = dimensionOp(a,2,scipvar.MUL);
                else %use n as dimension
                    if(nargin == 1)
                        n = 1;
                    end
                    c = dimensionOp(a,n,scipvar.MUL);
                end
            else
                c = a;
            end
        end

        %-- MISC UTILITIES --%
        
        %DISPLAY
        function display(a)
        %Overloaded Display method
           fprintf('\n------------------------------------------------\n');
           if(~isscalar(a))
               fprintf('\nSCIPVAR Object: %d x %d\n',size(a,1),size(a,2));
           else
               fprintf('\nScalar SCIPVAR Object\n\n');
           end
           strs = {'NUM','VAR','','MUL','DIV','ADD','SUB','SQUARE','SQRT','POW','EXPNT','LOG','SIN','COS','TAN','MIN','MAX','ABS','SIGN'};
           if(isscalar(a))
                %Print instruction tree
                for i = 1:2:length(a.ins)
                    fprintf('%-6s: %g\n',strs{a.ins(i)+1},a.ins(i+1));
                end                
           else
                [r,c] = size(a);
                for j = 1:c
                   for i = 1:r
                       fprintf('Eq(%d,%d) : %d instructions\n',i,j,length(a(i,j).ins)/2);
                   end
               end
           end
           fprintf('\n------------------------------------------------\n');
        end                
        
    end
    
    methods (Access = private)                
        
        %Basic Scalar and Vector Operations (+-*/^)
        function c = vectorOp(a,b,op)
            %Constants
            SCIP_NUM = 0;
            SCIP_VAR = 1;            
            %Ensure we have doubles (otherwise cast)
            if(isnumeric(a))
                if(~isa(a,'double'))
                    a = double(a);
                end
                if(issparse(a))
                    a = full(a);
                end
            end
            if(isnumeric(b))
                if(~isa(b,'double'))
                    b = double(b);
                end
                if(issparse(b))
                    b = full(b);
                end
            end
            %Now perform operation based on variable types presented
            switch [class(a),class(b)] 
                case 'scipvardouble'

                    %x1+2
                    if(isscalar(a) && isscalar(b))
                        c = a; %copy 
                        if(isempty(c.ins))
                            c.ins = [SCIP_VAR; c.indx];
                        end
                        c.ins = [c.ins;
                                 SCIP_NUM; b;
                                 op; NaN];

                    %[x1..xn]+2
                    elseif(isscalar(b))
                        c = a; %copy  
                        [r1,c1] = size(c);
                        for i = 1:r1
                            for j = 1:c1
                                if(isempty(c(i,j).ins))
                                    c(i,j).ins = [SCIP_VAR; c(i,j).indx];
                                end
                                c(i,j).ins = [c(i,j).ins;
                                              SCIP_NUM; b;
                                              op; NaN];
                            end
                        end
                    
                    %x1+[1..n]
                    elseif(isscalar(a))
                        %Create copies of current contents
                        [nrow,ncol] = size(b);
                        vins = cell(nrow,ncol);
                        vidx = cell(nrow,ncol);
                        if(isempty(a.ins))
                            a.ins = [SCIP_VAR; a.indx];
                        end
                        for i = 1:nrow
                            for j = 1:ncol
                                vins{i,j} = [a.ins; SCIP_NUM; b(i,j); op; NaN];
                                vidx{i,j} = a.indx;
                            end
                        end
                        %Create new object
                        c = scipvar(vins,vidx);

                    %[x1..xn]+[1..n]
                    else
                        %Check dimensions
                        [r1,c1] = size(a);
                        [r2,c2] = size(b);
                        if(r1 ~= r2)
                            error('Vector/Matrix row dimensions must agree');
                        elseif(c1 ~= c2)
                            error('Vector/Matrix column dimensions must agree');
                        end
                        c = a; %copy
                        for i = 1:r1
                            for j = 1:c1
                                if(isempty(c(i,j).ins))
                                    c(i,j).ins = [SCIP_VAR; c(i,j).indx];
                                end
                                c(i,j).ins = [c(i,j).ins;
                                              SCIP_NUM; b(i,j);
                                              op; NaN];
                            end
                        end                       
                    end
                    
                case 'doublescipvar'
                    %2+x1
                    if(isscalar(a) && isscalar(b))
                        c = b; %copy
                        if(isempty(c.ins))
                            c.ins = [SCIP_NUM; a
                                     SCIP_VAR; c.indx;
                                     op; NaN];                            
                        else
                            c.ins = [c.ins; 
                                     SCIP_NUM; a;
                                     op; 1]; % always flip (e.g. 2-2*x)
                        end
                        
                        
                    %2+[x1..xn]
                    elseif(isscalar(a))
                        c = b; %copy   
                        [r1,c1] = size(c);
                        for i = 1:r1
                            for j = 1:c1
                                if(isempty(c(i,j).ins))
                                    c(i,j).ins = [SCIP_NUM; a;
                                                  SCIP_VAR; c(i,j).indx;
                                                  op; NaN];
                                else
                                    c(i,j).ins = [c(i,j).ins;
                                                  SCIP_NUM; a;
                                                  op; 1]; %always flip
                                end
                            end
                        end
                        
                    %[1..n]+x1
                    elseif(isscalar(b))
                        %Create copies of current contents
                        [nrow,ncol] = size(a);
                        vins = cell(nrow,ncol);
                        vidx = cell(nrow,ncol);
                        if(isempty(b.ins))
                            b.ins = [SCIP_VAR; b.indx];
                        end
                        for i = 1:nrow
                            for j = 1:ncol
                                vins{i,j} = [b.ins; SCIP_NUM; a(i,j); op; 1];
                                vidx{i,j} = b.indx;
                            end
                        end
                        %Create new object
                        c = scipvar(vins,vidx);

                    %[1..n]+[x1..xn]
                    else
                        %Check dimensions
                        [r1,c1] = size(a);
                        [r2,c2] = size(b);
                        if(r1 ~= r2)
                            error('Vector/Matrix row dimensions must agree [r1 = %d, r2 = %d]',r1,r2);
                        elseif(c1 ~= c2)
                            error('Vector/Matrix column dimensions must agree [c1 = %d, c2 = %d]',c1,c2);
                        end
                        c = b; %copy
                        for i = 1:r1
                            for j = 1:c1
                                if(isempty(c(i,j).ins))
                                    c(i,j).ins = [SCIP_NUM; a(i,j);
                                                  SCIP_VAR; c(i,j).indx;
                                                  op; NaN];
                                else
                                    c(i,j).ins = [c(i,j).ins;
                                                  SCIP_NUM; a(i,j);
                                                  op; 1]; %always flip
                                end
                            end
                        end
                    end
                    
                case 'scipvarscipvar'
                    switch(op)
                        case scipvar.POW
                            %SCIP uses x^y = exp(y*log(x))
                            c = exp(b.*log(a));                            
                        otherwise   
                            %x1 + x2
                            if(isscalar(a) && isscalar(b))
                                c = a;
                                if(isempty(c.ins))
                                    c.ins = [SCIP_VAR; a.indx];
                                end
                                if(isempty(b.ins))
                                    b.ins = [SCIP_VAR; b.indx];
                                end
                                c.ins = [c.ins;
                                         b.ins;
                                         op; NaN];
                                     
                            %x1+[x2..xn]
                            elseif(isscalar(a))
                                c = b;
                                [r1,c1] = size(c);
                                if(isempty(a.ins))
                                    a.ins = [SCIP_VAR; a.indx];
                                end
                                for i = 1:r1
                                    for j = 1:c1
                                        if(isempty(c(i,j).ins))
                                            c(i,j).ins = [a.ins;
                                                          SCIP_VAR; c(i,j).indx;
                                                          op; NaN];
                                        else
                                            c(i,j).ins = [a.ins;
                                                          c(i,j).ins;                                                          
                                                          op; NaN]; 
                                        end
                                    end
                                end
                                
                            %[x1..n]+xk
                            elseif(isscalar(b))
                                c = a;
                                [r1,c1] = size(c);
                                if(isempty(b.ins))
                                    b.ins = [SCIP_VAR; b.indx];
                                end
                                for i = 1:r1
                                    for j = 1:c1
                                        if(isempty(c(i,j).ins))
                                            c(i,j).ins = [SCIP_VAR; c(i,j).indx;
                                                          b.ins;
                                                          op; NaN];
                                        else
                                            c(i,j).ins = [c(i,j).ins;
                                                          b.ins;
                                                          op; NaN]; 
                                        end
                                    end
                                end
                            else
                                %Check dimensions
                                [r1,c1] = size(a);
                                [r2,c2] = size(b);
                                if(r1 ~= r2)
                                    error('Vector/Matrix row dimensions must agree');
                                elseif(c1 ~= c2)
                                    error('Vector/Matrix column dimensions must agree');
                                end
                                c = a; %copy
                                for i = 1:r1
                                    for j = 1:c1
                                        if(isempty(c(i,j).ins))
                                            c(i,j).ins = [SCIP_VAR; c(i,j).indx];
                                        end
                                        if(isempty(b(i,j).ins))
                                            b(i,j).ins = [SCIP_VAR; b(i,j).indx];
                                        end
                                        c(i,j).ins = [c(i,j).ins;
                                                      b(i,j).ins;
                                                      op; NaN];                                        
                                    end
                                end
                            end
                    end
                    
                otherwise
                    switch(op)
                        case scipvar.ADD
                            error(['Can''t add (+) ',class(a),' and ',class(b)]);
                        case scipvar.SUB
                            error(['Can''t subtract (-) ',class(a),' and ',class(b)]);
                        case scipvar.DIV
                            error(['Can''t divide (/) ',class(a),' and ',class(b)]);
                        case scipvar.MUL
                            error(['Can''t multiply (*) ',class(a),' and ',class(b)]);
                        case scipvar.POW
                            error(['Can''t raise to a power (^) ',class(a),' and ',class(b)]);
                        otherwise
                            error(['Can''t (', op, ') ', class(a),' and ',class(b)]);
                    end
            end 
        end
        
        %Basic Vectorized Functions (log,log10,exp,abs)
        function c = vectorFcn(a,op)
           SCIP_VAR = 1;
           c = a;           
           if(isscalar(a))
                if(isempty(c.ins))
                    c.ins = [SCIP_VAR; a.indx];
                end
                c.ins = [c.ins; op; NaN];
           else
               [row,col] = size(c);
               for i = 1:row
                   for j = 1:col
                       if(isempty(c(i,j).ins))
                           c(i,j).ins = [SCIP_VAR; c(i,j).indx;
                                         op; NaN];
                       else
                           c(i,j).ins = [c(i,j).ins;
                                         op; NaN];
                       end
                   end
               end
           end            
        end 
        
        %Dimension Wise Operations (sum, prod)
        function b = dimensionOp(a,n,op)            
            [r,c] = size(a);
            if(n == 1) %column products / sums
                b = a(1,1);
                if(c == 1) %result is a scalar
                    for i = 2:r
                        switch(op)
                            case scipvar.ADD, b = b + a(i,1);
                            case scipvar.MUL, b = b * a(i,1);
                        end
                    end
                else %result is a row vector
                    b = a(1,:);
                    for i = 2:r
                        for j = 1:c
                            switch(op)
                                case scipvar.ADD, b(j) = b(j) + a(i,j);
                                case scipvar.MUL, b(j) = b(j) * a(i,j);
                            end
                        end
                    end
                end               
            elseif(n==2) %row products / sums
                b = a(1,1);
                if(r == 1) %result is a scalar
                    for i = 2:c
                        switch(op)
                            case scipvar.ADD, b = b + a(1,i);
                            case scipvar.MUL, b = b * a(1,i);
                        end
                    end
                else %result is a column vector
                    b = a(:,1);
                    for i = 1:r
                        for j = 2:c
                            switch(op)
                                case scipvar.ADD, b(i) = b(i) + a(i,j);
                                case scipvar.MUL, b(i) = b(i) * a(i,j);
                            end
                        end
                    end
                end
            else
                error('Only 2D operations are implemented');
            end           
        end
        
        %Matrix Operations
        function c = matrixOp(a,b,op)            
            %Check dimensions
            [r1,c1] = size(a);
            [r2,c2] = size(b);               
            if(c1 ~= r2)
                error('Matrix / vector sizes don''t correspond - trying to multiply %d x %d by %d x %d',r1,c1,r2,c2);
            end            
            %Create new output object, then fill it
            c = scipvar(r1,c2);
            for i = 1:r1
                for j = 1:c2
                    c(i,j) = dot(a(i,:),b(:,j)');
                end
            end
        end
    end    
end

