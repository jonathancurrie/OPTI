function p = sdpWrite(prob,filename,type)
%sdpWrite  Write a Semidefinite Programming File from Matlab Values
%
% sdpWrite(prob,filename,type) writes the OPTI SDP problem to an SDPA or 
% SeDuMi file. If you do not specify a file extension it will default to 
% the type specified.
%
%   Available File Types:
%       - SDPA-S        [.dat-s]        (Sparse Format)
%       - SDPA          [.dat]          (Dense Format)
%       - SeDuMi        [.mat]          
%
% You may specify a full path to the file, or if you specify a filename
% only, it will be written to the current directory.

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
            error('Unknown file extension: %s. Please specify the file type using a second arugment');
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

%Ensure constraints in general format with only inequalities
if(~isempty(prob.ru))
    [A,b,Aeq,beq] = row2gen(prob.A,prob.rl,prob.ru);
    [prob.A,prob.b] = inequal(A,b,Aeq,beq);
end

%Call read function based on file type
switch(filetype)
    case {'SDPA','SDPA-S'}
        %Check if we are in SeDuMi format
        if(isstruct(prob.sdcone))
            prob = sedumi2opti(prob);
        end
        if(isstruct(prob.int))
            prob.int = prob.int.str;
        end
        p = optiWriteSDPA(filename,prob.f,prob.A,prob.b,prob.lb,prob.ub,prob.sdcone,prob.int,strcmpi(filetype,'SDPA'));
    case 'SEDUMI'
        [At,b,c,K] = opti2sedumi(prob); %#ok<NASGU,ASGLU>
        save(filename,'At','b','c','K');
end


