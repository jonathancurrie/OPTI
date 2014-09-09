function coinWrite(prob,filename,type)
%COINWRITE  Write a Mathematical File From Matlab Values
%
% coinWrite(prob,filename,type) writes the problem prob to the file specified 
% by filename and problem type specified by type. If you do not specify a 
% file extension it will default to the type specified.
%
%   Available File Types:
%       - MPS  [.mps]
%       - QPS  [.qps]
%       - LP   [.lp]  (Appears to be CPLEX LP version - Simplified)
%
% You may specify a full path to the file, or if you specify a filename
% only, it will be written to the current MATLAB directory.
%
% The routines underneath use COIN-OR utilities for File IO. See 
% attached EPL License.

%   Copyright (C) 2011 Jonathan Currie (I2C2)

%Quick Checks
if(~exist('type','var') || isempty(type))
    %See if contained within file name
    dots = strfind(filename,'.');
    if(isempty(dots))
        error('You must supply the file type to this function');
    else %hopefully got a good ext, will check below
        type = filename(dots(end)+1:end);
    end
end
if(~ischar(filename))
    error('Filename must be a char array!');
end

%Problem Checks
if(~isempty(prob.fun) || ~isempty(prob.nlcon))
    error('You cannot write Nonlinear problems using this interface');
end
if(isa(prob.f,'function_handle') || isa(prob.H,'function_handle'))
    error('f and H must be vectors / matrices!');
end
if(~isempty(prob.Q))
    error('You cannot write quadratic constraints using this interface');
end
if(prob.sense == -1)
    error('You cannot explicity write maximization problems using this interface'); %fix up later
end

%Determine file type
switch(lower(type))
    case 'mps'
        filetype = 'mps';
    case 'qps'
        filetype = 'qps';
    case 'lp'
        filetype = 'lp';
    otherwise
        error('Unknown file type %s',type);
end

%Check problem type vs file type
if(isfield(prob,'type') && ~isempty(prob.type))
    switch(lower(prob.type))
        case {'lp','bilp','milp'}
            if(strcmpi(filetype,'qps'))
                error('LP/BILP/MILPs should be written to MPS or LP files');
            end
        case {'qp','miqp'}
            if(~strcmpi(filetype,'qps'))
                error('QP/MIQPs should be written to QPS files');
            end
        otherwise
            error('Problem type %s is not supported by this interface',upper(prob.type));
    end
end         

% Make sure we have file extension
if(isempty(strfind(filename,['.' filetype])))
    filename = [filename '.' filetype];
end 

% If filename includes absolute path, we can skip which
if(~isempty(strfind(filename,':')))
    p = filename;
else %Locate the full path to the file
    p = [cd filesep filename];
end

%Convert to coin problem
if(~isempty(prob.b) || ~isempty(prob.beq))
    [prob.A,prob.rl,prob.ru] = gen2row(prob.A,prob.b,prob.Aeq,prob.beq);
end
if(isempty(prob.A))
    prob.A = spalloc(0,length(prob.f),0); %error thrown internally if not ok
end

%Ensure sparsity
prob.A = sparse(prob.A);
prob.H = sparse(prob.H);

%Ensure Bounds
if(isempty(prob.lb))
    prob.lb = -Inf(size(prob.f));
end
if(isempty(prob.ub))
    prob.ub = Inf(size(prob.f));
end

%Setup Integer Vars
if(isstruct(prob.int))
    prob.int = prob.int.str;
end 
if(isempty(prob.int))
    prob.int = zeros(size(prob.f));
elseif(ischar(prob.int)) %char string
    prob.int = lower(prob.int);
    i = zeros(size(prob.f));
    int = prob.int == 'i';
    bin = prob.int == 'b';
    i(int) = 1;
    i(bin) = 1;
    prob.lb(bin) = 0;
    prob.ub(bin) = 1;
    prob.int = i;
elseif(length(prob.int) ~= length(prob.f) || any(prob.int > 1)) %find indicies (might be wrong though!)
    i = zeros(size(prob.f));
    i(prob.int) = 1;
    prob.int = i;
end
prob.int = int8(prob.int);

%Setup SOS
if(~isempty(prob.sos))
    [r,c] = size(prob.sos.type);
    if(r > c)
        prob.sos_type = str2num(prob.sos.type); %#ok<ST2NM>
    else
        prob.sos_type = str2num(prob.sos.type'); %#ok<ST2NM>
    end
    if(length(prob.sos.type) == 1)
        prob.sos_index = {prob.sos.index};
        prob.sos_weight = {prob.sos.weight};
    else
        prob.sos_index = prob.sos.index;
        prob.sos_weight = prob.sos.weight;
    end
else
    prob.sos_type = [];
    prob.sos_index = [];
    prob.sos_weight = [];
end
%Setup Objc
if(~isempty(prob.objbias) && prob.objbias == 0)
    prob.objbias = 0;
end

%Call coin-or mex
coinW(prob,p,filetype);
