function [prob,p] = sdpRead(filename,type,print)
%sdpRead  Read a Semidefinite Programming File and convert to Matlab Values
%
% prob = sdpRead(filename,type) reads the file specified by filename
% and problem type specified by type and converts it to an optiprob. If you
% do not specify a file extension it will default to the type specified.
% The returned structure is solver independent so the user can manually 
% extract matrices as required or supply it directly to opti().
%
% prob = sdpRead(filename,type,print) specifies to display information
% and error messages in MATLAB when print > 0. The default is 0.
%
%   Available File Types:
%       - SDPA-S        [.dat-s]        (Sparse Format)
%       - SDPA          [.dat]          (Dense Format)
%       - SeDuMi        [.mat]          
%
% You may specify a full path to the file, or if you specify a filename
% only, it must exist on the MATLAB path.

%   Copyright (C) 2013 Jonathan Currie (IPL)

%Quick Checks
if(~exist('type','var') || isempty(type))
    %See if contained within file name
    dots = strfind(filename,'.');
    if(isempty(dots))
        error('You must supply the file type to this function');
    else %hopefully got a good ext, will check below
        ext = filename(dots(end)+1:end);
    end
else
    ext = [];
end
if(~exist('print','var') || isempty(print))
    print = 0;
end
if(~ischar(filename))
    error('Filename must be a char array!');
end

%If we have a file extension and no type, use extension to determine file
%type
if(~isempty(ext))
    switch(lower(ext))
        case 'dat-s'
            filetype = 'SDPA-S';
        case 'dat'
            filetype = 'SDPA';
        case 'mat'
            filetype = 'SEDUMI';
        otherwise
            error('Unknown file extension: %s. Please specify the file type using a second arugment',ext);
    end
else
    %Determine file type from user specification
    switch(lower(type))
        case 'sdpa-s'
            filetype = 'SDPA-S';
            ext = 'dat-s';
        case 'sdpa'
            filetype = 'SDPA';
            ext = 'dat';
        case 'sedumi'
            filetype = 'SEDUMI';
            ext = 'mat';
        otherwise
            error('Unknown file type %s',type);
    end
end

% Make sure we have file extension
if(isempty(strfind(lower(filename),['.' ext])))
    filename = [filename '.' ext];
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

%Call read function based on file type
switch(filetype)
    case {'SDPA','SDPA-S'}
        [f,A,b,sdcone,xtype] = optiReadSDPA(filename,strcmpi(filetype,'SDPA'),print);
        prob = optiprob('f',f,'ineq',A,b,'sdcone',sdcone,'xtype',xtype);
    case 'SEDUMI'
        load(filename);
        if(~exist('At','var')), error('Supplied SeDuMi file does not contain the variable At'); end
        if(~exist('b','var')),  error('Supplied SeDuMi file does not contain the variable b'); end
        if(~exist('c','var')),  error('Supplied SeDuMi file does not contain the variable c'); end
        if(~exist('K','var')),  error('Supplied SeDuMi file does not contain the variable K'); end
        prob = optiprob('sedumi',At,full(b),full(c),K);
end


