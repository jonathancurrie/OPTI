function opti_package
% OPTI_PACKAGE  Generates GitHub releases .zip files

% Get Solver Names (ones we want to upload)
s = lower(optiSolver('all'));
removeIdx = strcmp(s,'baron') | strcmp(s,'cplex') | strcmp(s,'gmatlab') | strcmp(s, 'matlab') | strcmp(s,'mosek') | strcmp(s,'scip') | strcmp(s,'sedumi');
solvers = s(~removeIdx);

% Utility Names
utilities = {'asl', 'coinR', 'coinW', 'mklJac', 'rmathlib'};

% Build zip file list
files = [solvers utilities];
for i = 1:length(files)
    files{i} = which([files{i} '.' mexext]);
end

% Zip the files
fprintf('Compressing MEX Files...');
name = ['optiMEXFiles_' mexext '_' sprintf('%.0f_%.0f',optiver,mod(optiver*100,100))];
zip(name,files);
rehash;
rehash;
fprintf('Done!\n');

% Update the contents file
UpdateContentsFile('OPTI',optiver,'Utilities/opti');

% Copy latest SCIP to download folder
copyfile(which(['scip.' mexext]), '../../../OPTI Academic Solvers');

function UpdateContentsFile(name,tbxver,contentsFile)
%Update contents file description line

if(~isempty(strfind(contentsFile,'Contents.m')))
    p = contentsFile;
elseif(contentsFile(end) == filesep)
    p = [contentsFile 'Contents.m'];
else
    p = [contentsFile filesep 'Contents.m'];
end
s = fileread(p);
%Find initial date
idx = strfind(s,'(C)');
idxe = strfind(s(idx:end),'-');
stDate = strtrim(s(idx+3:idx+idxe-2));
mver = ver('matlab');
%Regenerate Contents File
fid = fopen(p,'w');
if(fid)
    try
        fprintf(fid,'%% %s Toolbox\n',upper(name));
        fprintf(fid,'%% Version %.2f %s %s\n',tbxver,mver.Release,datestr(now,1));
        fprintf(fid,'%% Copyright (C) %s-%s Jonathan Currie (Control Engineering)\n',stDate,datestr(now,10));
        fprintf(fid,'%% License: https://controlengineering.co.nz/Wikis/OPTI/index.php/DL/License\n');
    catch ME
        fclose(fid);
        rethrow(ME);
    end
end
fclose(fid);