function [source,header] = VS_BuildFileList(directory,exFolder,mode)
% Create a file list from a selected source directory.

% Horribly memory inefficient code! Will think about rewriting if this
% actually becomes useful to anyone.

source = {}; j = 1;
header = {}; k = 1;

if(nargin < 2 || isempty(exFolder))
    exFolder = {};
end
if(nargin < 3 || isempty(mode))
    mode = 1; % c/c++
end

%For each path passed, read in files
for m = 1:size(directory,1)
    %If multiple paths, read each
    if(iscell(directory))
        d = directory{m};
    else
        d = directory;
    end

    %Get all folders under the directory
    p = genpath(d);

    if(isempty(p))
        error(['Cannot find directory: ' d]);
    end

    %Split into directories
    ps = regexp(p,';','split');
    
    %Remove unwanted folders
    if(~isempty(exFolder))
        len = length(d)+2;
        for i = 1:length(ps)
            if(length(ps{i}) > len-2)
                if(ismember(ps{i}(len:end),exFolder))
                    ps{i} = [];
                end
            end
        end
        %Try again, this time just searching end folder
        for i = 1:length(ps)
            ind = strfind(ps{i},filesep);
            if(~isempty(ind))
                pa = ps{i}(ind(end):end);
                for j = 1:length(exFolder)
                    if(strcmp(pa,[filesep exFolder{j}]))
                        ps{i} = [];
                        break;
                    end
                end
            end
        end
    end

    for i = 1:length(ps)
        %Read Files in Current Path
        if(~isempty(ps{i}))
            [src,hdr,edir] = readFiles(ps{i},mode);
            if(~isempty(src) || ~isempty(edir))
                source{j,1} = ps{i};
                source{j,2} = src;
                source{j,3} = m;
                j = j + 1;
            end
            if(~isempty(hdr) || ~isempty(edir))
                header{k,1} = ps{i};
                header{k,2} = hdr;
                header{k,3} = m;
                k = k + 1;
            end
        end
    end
end


function [src,hdr,edir] = readFiles(p,mode)
%Build a cell array containing source files

src = {}; j = 1;
hdr = {}; k = 1;
edir = [];
f = dir(p);

if(isempty(f)); return; end;

for i = 1:length(f)
    if(~f(i).isdir)  
        n = f(i).name;
        switch(mode)
            case 1
                %C/C++ Files
                if(length(n) > 3)
                    if(strcmpi(n(end-1:end),'.c') || strcmpi(n(end-3:end),'.cpp') || strcmpi(n(end-2:end),'.cc'))
                        src{j,1} = n; %#ok<*AGROW>
                        j = j + 1;
                    end
                elseif(length(n) > 2)
                    if(strcmpi(n(end-1:end),'.c'))
                        src{j,1} = n;
                        j = j + 1;
                    end
                end
                %Header Files
                if(length(n) > 3)
                    if(strcmpi(n(end-1:end),'.h') || strcmpi(n(end-3:end),'.hpp'))
                        hdr{k,1} = n;
                        k = k + 1;
                    end
                elseif(length(n) > 2)
                    if(strcmpi(n(end-1:end),'.h'))
                        hdr{k,1} = n;
                        k = k + 1;
                    end
                end
            case 0
                %Fortran Files
                if(strcmpi(n(end-1:end),'.f'))
                    src{j,1} = n;
                    j = j + 1;
                end
        end
    end
end

%Directory has no files, but still required for subdirectories
if(isempty(src) && isempty(hdr) && (size(f,1) > 2))
    edir = 1;
end



