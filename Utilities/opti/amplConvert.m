function filename = amplConvert(filename,extraArgs,amplPath)
%AMPLCONVERT  Read an AMPL .MOD File and convert it to .NL using AMPL
%
% prob = amplRead(filename) will read the .mod file and attempt to find
% ampl.exe to convert it to a .nl. The AMPL executable must be on the
% MATLAB path. The resulting .nl file will be in the same directory as the
% original .mod.
%
% prob = amplRead(filename,extraArgs) allows you to pass extra
% arguments to AMPL when converting your model file, such as options or
% data. Files passed must include the file extension, and full path if not
% on the MATLAB path. Pass multiple arguments as a cell array.
%
% prob = amplRead(filename,extraArgs,amplPath) reads the .MOD (AMPL Model) 
% specified by filename, then converts it using your AMPL executable 
% (specified by amplPath) to a .NL file. filename must include the full 
% path if not on the MATLAB path.
% 
% You must have a licensed version of AMPL (or the student edition) present
% on your computer. For more information see www.ampl.com.

%   Copyright (C) 2013 Jonathan Currie (I2C2)

% Attempt to find AMPL exectuable
if(nargin < 3 || isempty(amplPath))
    amplPath = which('ampl.exe');
    if(isempty(amplPath))
        error('Cannot find the AMPL executable (ampl.exe) on the MATLAB path. Please supply the AMPL path as the third argument to this function (or amplRead).');
    end
end
% See if AMPL path contains the executable
if(~ischar(amplPath))
    error('amplPath must be a char array!');
end
if(isempty(strfind(amplPath,'.exe')))
    %Assume called ampl.exe and tag on
    if(amplPath(end) == filesep)
        amplPath = [amplPath 'ampl.exe'];
    else
        amplPath = [amplPath filesep 'ampl.exe'];
    end
end
%Check finally we have the entire AMPL string
if(exist(amplPath,'file') == 0)
    error('Cannot find AMPL executable. Checked the following path:\n%s',amplPath);
end

%Check we don't have a .nl
if(~isempty(strfind(filename,'.nl')))
    error('This function is for converting .mod files to .nl only');
end
% Make sure we have file extension
if(isempty(strfind(filename,'.mod')))
    filename = [filename '.mod'];
end 
% If filename includes absolute path, we can skip which
if(~isempty(strfind(filename,':')))
    mp = filename;
    %Check it exists
    if(~exist(filename,'file'))
        error('Cannot locate file %s!',filename);
    end
else %Locate the full path to the file (must be on matlab path)
    mp = which(filename);
    %Check it exists
    if(isempty(mp))
        error('Cannot locate file %s!',filename);
    end
end
%Get last full stop
e = strfind(mp,'.');
if(length(e) > 1)
    e = e(end); %assume is the end one?
end
fn = mp(1:e-1);

%Call AMPL, processing optional arguments if present
if(nargin < 2 || isempty(extraArgs))
    [status,result] = system(['"' amplPath '" -og"' fn '" "' mp '"']);
else
    %Build extra arg string
    post = [];
    if(~iscell(extraArgs))
        extraArgs = {extraArgs};
    end
    for i = 1:length(extraArgs)
        str = extraArgs{i};
        if(~ischar(str))
            error('Extra argument %d was not a string!',i);
        end
        %Check if we have a file (assume it includes the extension, so
        %check for a .
        if(~isempty(strfind(str,'.')))
            % If filename includes absolute path, we can skip which
            if(~isempty(strfind(str,':')))
                p = str;
                %Check it exists
                if(~exist(str,'file'))
                    error('Cannot locate extra arg file %s!',str);
                end
            else %Locate the full path to the file (must be on matlab path)
                p = which(str);
                %Check it exists
                if(isempty(p))
                    error('Cannot locate extra arg file %s!',str);
                end
            end
        else
            p = str; %take as is (could be a parameter)
        end
        post = [post ' "' p '" ']; %#ok<AGROW>
    end
    [status,result] = system(['"' amplPath '" -og"' fn '" "' mp '" ' post]);
end

if(status ~= 0)
    error('There was an error processing your AMPL file. Error:\n\n%s',result);
end

%Return filename appended with .nl
filename = [fn '.nl'];
