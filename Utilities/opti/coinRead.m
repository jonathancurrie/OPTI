function [prob,p] = coinRead(filename,type,print)
%coinRead  Read a Mathematical File and convert to Matlab Values
%
% prob = coinRead(filename,type) reads the file specified by filename
% and problem type specified by type and converts it to an optiprob. If you
% do not specify a file extension it will default to the type specified.
% The returned structure is solver independent so the user can manually 
% extract matrices as required or supply it directly to opti().
%
% prob = coinRead(filename,type,print) specifies to display information
% and error messages in MATLAB when print > 0. The default is 0.
%
%   Available File Types:
%       - MPS  [.mps]
%       - QPS  [.qps]
%       - LP   [.lp]  (Appears to be CPLEX LP version - Simplified)
%       - GMPL [.mod] (Subset of AMPL, implemented via GLPK)
%       - GAMS [.gms] (Only first COIN version... Will probably not work)
%
% You may specify a full path to the file, or if you specify a filename
% only, it must exist on the MATLAB path.
%
% The routines underneath use COIN-OR & GLPK utilities for File IO. See 
% attached EPL License for COIN-OR and GPL for GLPK.

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
if(~exist('print','var') || isempty(print))
    print = 0;
end
if(~ischar(filename))
    error('Filename must be a char array!');
end

%Determine file type
switch(lower(type))
    case 'mps'
        filetype = 'mps';
    case 'qps'
        filetype = 'qps';
    case 'lp'
        filetype = 'lp';
    case {'gmpl','mod'}
        filetype = 'mod';
    case {'gams','gms'}
        filetype = 'gms';
    otherwise
        error('Unknown file type %s',type);
end

%Check if extension is specified
if(isempty(strfind(filename,'.')))
    filename = [filename '.' filetype];
end

% If filename includes absolute path, we can skip which
if(~isempty(strfind(filename,':')))
    p = filename;
    %Check it exists
    if(~exist(filename,'file'))
        error('Cannot locate file %s!',filename);
    end
else %Locate the full path to the file (must be on matlab path)
    p = which(filename);
    %Check it exists
    if(isempty(p))
        error('Cannot locate file %s!',filename);
    end
end

%Call coin-or mex
[f,A,rl,ru,lb,ub,ivars,H,name,sostype,sosind,soswt,objc] = coinR(p,filetype,print);

%No guarantee matrices have been created with correct order (T.Kelman, 6/9/12), double transpose to ensure order)
H = H'';
A = A'';

%Remove H if all zero
if(~nnz(H))
    H = [];
%Assume we need to make symmetric if e.g. [8 0; 2 10]
elseif(nnz(triu(H,1)) == 0 && ~(nnz(tril(H,-1)) == 0))
    optiwarn('OPTI:COINR_QUADOBJ','Assumed QOBJ has been specified as lower triangular only, forcing to symmetric');
    H = H + tril(H,-1)';
%Check for upper tri only
elseif(nnz(tril(H,-1)) == 0 && ~(nnz(triu(H,1)) == 0))
    optiwarn('OPTI:COINR_QUADOBJ','Assumed QOBJ has been specified as upper triangular only, forcing to symmetric');
    H = H + triu(H,1)';
end
    
%Build Integer String
int = repmat('C',size(f')); ivars = logical(ivars);
int(ivars) = 'I'; %assume just integer variables (no binary - even with MPS)
%Build SOS
if(~isempty(sostype))
    sostype = num2str(sostype);
    if(length(sostype) == 1)
        sosind = sosind{1}; %remove cells
        soswt = soswt{1};
    end
else
    sostype = []; sosind = []; soswt = [];
end


%Build return problem
prob = optiprob('H',H,'f',f,'objbias',objc,'lin',A,rl,ru,'bounds',lb,ub,'int',int,'name',name,'sos',sostype,sosind,soswt);
