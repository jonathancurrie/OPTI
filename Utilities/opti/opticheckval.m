classdef opticheckval
    %OPTICHECKVAL Contains Static Methods for Checking Option Values
    
    methods(Static)
        
        function err = checkScalarSet(value,field,set)
            %Scalar in set
            err = [];
            if(~isscalar(value) || ~isnumeric(value) || ~isa(value,'double') || ~any(value==set))
                str = sprintf('\n');
                for i = 1:length(set)
                    str = sprintf('%s %g\n',str,set(i));
                end
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should a double scalar from one of the following supported values:\n%s',field,str);
            end
        end
        
        function err = checkScalar01(value,field)
            %Scalar 0 or 1
            err = [];
            if(~isscalar(value) || ~isnumeric(value) || (value ~= 0 && value ~=1) || ~isa(value,'double'))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be an integer scalar with the value 0 or 1',field);
            end
        end
        
        function err = checkVector01(value,field)
            %Vector 0 or 1
            err =[];
            if(~isnumeric(value) || any((value == 0 | value == 1) ~= 1) || ~isa(value,'double'))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a vector with the values 0 or 1',field);
            end
        end
        
        function err = checkScalarGrtZ(value,field)
            %Scalar double > Zero
            err = [];
            if(~isscalar(value) || ~isnumeric(value) || value <= 0 || ~isa(value,'double'))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a double scalar greater than 0',field);
            end
        end
        
        function err = checkScalarIntGrtZ(value,field)
            %Scalar integer > zero
            err = [];
            if(~isscalar(value) || ~isnumeric(value) || value < 1 || (round(value) ~= value) || ~isa(value,'double'))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be an integer scalar greater than 0',field);
            end
        end
        
        function err = checkScalarNonNeg(value,field)
            %Scalar double >= Zero
            err = [];
            if(~isscalar(value) || ~isnumeric(value) || value < 0 || ~isa(value,'double'))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a scalar greater or equal to 0',field);
            end            
        end   
        
        function err = checkScalarIntNonNeg(value,field)
            %Scalar integer >= zero
            err = [];
            if(~isscalar(value) || ~isnumeric(value) || value < 0 || (round(value) ~= value) || ~isa(value,'double'))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be an integer scalar greater than 0',field);
            end
        end
        
        function err = checkScalarDbl(value,field)
            %Scalar double
            err = [];
            if(~isscalar(value) || ~isnumeric(value) || ~isreal(value) || ~isa(value,'double'))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a scalar double',field);
            end
        end
        
        function err = checkScalarBoundLELE(value,field,lb,ub)
            %Scalar double lb <= value <= ub
            err = [];
            if(~isscalar(value) || ~isnumeric(value) || ~isreal(value) || value < lb || value > ub || ~isa(value,'double'))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a scalar double %g <= value <= %g',field,lb,ub); 
            end
        end
        
        function err = checkScalarBoundLLE(value,field,lb,ub)
            %Scalar double lb < value <= ub
            err = [];
            if(~isscalar(value) || ~isnumeric(value) || ~isreal(value) || value <= lb || value > ub || ~isa(value,'double'))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a scalar double %g < value <= %g',field,lb,ub); 
            end
        end
        
        function err = checkScalarBoundLEL(value,field,lb,ub)
            %Scalar double lb <= value < ub
            err = [];
            if(~isscalar(value) || ~isnumeric(value) || ~isreal(value) || value < lb || value >= ub || ~isa(value,'double'))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a scalar double %g <= value < %g',field,lb,ub); 
            end
        end
        
        function err = checkScalarBoundLL(value,field,lb,ub)
            %Scalar double lb < value < ub
            err = [];
            if(~isscalar(value) || ~isnumeric(value) || ~isreal(value) || value <= lb || value >= ub || ~isa(value,'double'))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a scalar double %g < value < %g',field,lb,ub); 
            end
        end
        
        function err = checkScalarIntBoundLEL(value,field,lb,ub)
            %Scalar integer lb <= value < ub
            err = [];
            if(~isscalar(value) || ~isnumeric(value) || ~isreal(value) || (round(value) ~= value) || value < lb || value >= ub || ~isa(value,'double'))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a scalar integer %g <= value < %g',field,lb,ub); 
            end
        end
        
        function err = checkScalarIntBoundLELE(value,field,lb,ub)
            %Scalar integer lb <= value <= ub
            err = [];
            if(~isscalar(value) || ~isnumeric(value) || ~isreal(value) || (round(value) ~= value) || value < lb || value > ub || ~isa(value,'double'))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a scalar integer %g <= value <= %g',field,lb,ub); 
            end
        end
        
        function err = checkScalarIntM1GrtZ(value,field)
            %Scalar integer value = -1 or > 0
            err = [];
            if(~isscalar(value) || ~isnumeric(value) || ~isreal(value) || (round(value) ~= value) || ~isa(value,'double') || (value ~= -1 && value < 1))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a scalar integer = -1 OR > 0',field); 
            end
        end
        
        function err = checkVectorIntGrtZ(value,field)
            %Vector of integers > 0
            err = [];
            if(~isnumeric(value) || any(value < 1) || any(round(value) ~= value) || ~any(isa(value,'double')))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be an integer vector with all elements greater than 0',field);
            end
        end
        
        function err = checkDualCol(value,field)
            %Double matrix with two columns
            err = [];
            if(~isnumeric(value) || ~isreal(value) || size(value,2) ~= 2 || ~isa(value,'double'))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a double matrix with two columns',field);
            end
        end
        
        function err = checkDblVec(value,field)
            %Double vector (not sparse)
            err = [];
            if(~isnumeric(value) || ~isreal(value) || ~isa(value,'double') || (size(value,1) > 1 && size(value,2) > 1) || issparse(value))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a real double dense vector (1D)',field);
            end            
        end
        
        function err = checkDblMat(value,field)
            %Double Matrix or Vector or Scalar
            err = [];
            if(~isnumeric(value) || ~isreal(value) || ~isa(value,'double'))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a real double matrix or vector (2D or 1D)',field);
            end            
        end
        
        function err = checkDblVecOrCell(value,field)
            %Double Vector or Cell Array of Double Vectors
            err = [];
            if(~iscell(value)), value = {value}; end
            for i = 1:length(value)
                val = value{i};
                if(~isnumeric(val) || ~isreal(val) || ~isa(val,'double') || (size(val,1) > 1 && size(val,2) > 1) || issparse(val))
                    err = MException('OPTI:SetFieldError','Parameter ''%s'' (cell %d) should be a real double dense vector or scalar (1D)',field,i);
                end
            end
        end
        
        function err = checkDblVecOrLogVec(value,field)
            %Double Vector or Logical Vectors
            err = [];
            if(islogical(value))
                if(size(value,1) > 1 && size(value,2) > 1)
                    err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a logical vector OR real double dense vector (1D)',field);
                end
            else
                if(~isnumeric(value) || ~isreal(value) || ~isa(value,'double') || (size(value,1) > 1 && size(value,2) > 1) || issparse(value))
                    err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a real double dense vector (1D) or logical vector',field);
                end 
            end
        end
        
        function err = checkDblMatOrCell(value,field)
            %Double Matrix or Cell Array of Double Matrices
            err = [];
            if(~iscell(value)), value = {value}; end
            for i = 1:length(value)
                val = value{i};
                if(~isnumeric(val) || ~isreal(val) || ~isa(val,'double'))
                    err = MException('OPTI:SetFieldError','Parameter ''%s'' (cell %d) should be a real double matrix or vector (2D or 1D)',field,i);
                end
            end
        end
        
        function err = checkChar(value,field)
            %Char array
            err = [];
            if(~ischar(value))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a char array (string)',field);
            end  
        end
        
        function err = checkNLOPTAlg(value,field)
            %Char array and valid NLOPT algorithm
            err = [];
            if(~ischar(value))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a char array (string)',field);
            else
                nloptSolver(value); %will throw an error if not found
            end
        end
        
        function err = checkOptiSolver(value,field)
            %Char array and valid OPTI solver
            err = [];
            if(~ischar(value))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a char array (string)',field);
            else
                rub = checkSolver(value,1); %#ok<NASGU> 
            end
        end
        
        function err = checkYesNo(value,field)
            %Char array and yes or no
            err = opticheckval.checkValidString(value,field,{'yes','no'});
        end
        
        function err = checkYesCrashNo(value,field)
            %Char array and no, or yes will crash
            err = [];
            if(~ischar(value) || ~strcmpi(value,'no'))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' must be set as no, otherwise OPTI may crash (sorry!).\n',field);
            end
        end
        
        function err = checkEnDis(value,field)
            %Char array and enable or disable
            err = opticheckval.checkValidString(value,field,{'enable','disable'});
        end
        
        function err = checkOnOff(value,field)
            %Char array and on or off
            err = opticheckval.checkValidString(value,field,{'on','off'});
        end
        
        function err = checkValidString(value,field,valid)
            %Char array and valid string
            err = [];
            if(~ischar(value) || all(strcmpi(value,valid)==0))
                str = sprintf('\n');
                for i = 1:length(valid)
                    str = sprintf('%s ''%s''\n',str,valid{i});
                end
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should a string from the following supported options:\n%s',field,str);
            end
        end
        
        function err = checkValidStringArr(value,field,valid)
            %Char array and of valid characters
            err = [];
            if(~ischar(value))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a char array (string)',field);
            else
                if(~iscell(valid)), valid = {valid}; end
                for i = 1:length(value)
                   if(~any(ismember(valid,value(i))))
                       str = sprintf('\n');
                        for j = 1:length(valid)
                            str = sprintf('%s ''%s''\n',str,valid{j});
                        end
                        err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a char array (string) with the following characters only:\n%s',field,str);
                        return;
                    end
                end
            end
        end
        
        function err = checkStruct(value,field)
            %Structure
            err = [];
            if(~isstruct(value))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a structure\n',field);
            end
        end
        
        function err = checkStructFields(value,field,reqfields)
            %Structure with required fields
            err = [];
            if(~isstruct(value))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a structure\n',field);
            else
                flds = fieldnames(value);
                if(~iscell(reqfields)), reqfields = {reqfields}; end
                for i = 1:length(reqfields)
                    if(~any(ismember(flds,reqfields{i})))
                        str = sprintf('\n');
                        for j = 1:length(reqfields)
                            str = sprintf('%s ''%s''\n',str,reqfields{j});
                        end
                        err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a structure with the following fields:\n%s',field,str);
                        return;
                    end
                end
            end
        end

        function err = checkFunHandle(value,field)
            %Function Handle
            err = [];
            if(~isa(value,'function_handle'))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a function handle.\n\nFor existing m-file functions, pass as @function_name (not a string, using the ''@'' operator).\n',field);
            end
        end
        
        function err = checkNumorFHndl(value,field)
            %Function handle or numeric matrix
            err = [];
            if(isnumeric(value))
                if(~isreal(value) || ~isa(value,'double'))
                    err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a real, double, matrix/vector or function handle\n',field);
                end
            elseif(~isa(value,'function_handle'))
                err = MException('OPTI:SetFieldError','Parameter ''%s'' should be a function handle OR real, double, matrix\n',field);
            end
        end
        
        
        function err = relErrorCheck(v1,v2,str,tol)
            %Relative Error Check between two scalars (v1 calc, v2 true)
            if(nargin < 4 || isempty(tol)), tol = 1e-6; end
            if(nargin < 3), str = []; end
            
            if(v2 < eps(class(v2)))
                err = 0;
            elseif(v1==0 || v2==0)
                err = abs(v1-v2);
            else
                err = abs(v1-v2)/v1;
            end
            if(err > tol)
                throwAsCaller(MException('OPTI:RelativeErrorCheck','Error when comparing parameter ''%s'' - Relative Error %g%%\n',str,err*100));
            end
        end
    end
    
end

