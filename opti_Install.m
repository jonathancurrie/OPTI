function opti_Install
%% Installation File for OPTI

% In order to run this tool, please run this file to setup the required
% directories. You MUST be in the current directory of this file!

%   Copyright (C) 2016 Jonathan Currie (IPL)

cpath = cd;
try
    cd('Utilities/opti');
catch %#ok<CTCH>
    error('You don''t appear to be in the OPTI Toolbox directory');
end
%Get current versions    
localVer = optiver();

fprintf('\n------------------------------------------------\n')
fprintf(['  INSTALLING OPTI TOOLBOX ver ' sprintf('%1.2f',localVer) '\n'])

cd(cpath);
% Check ML ver
matlabVerCheck();

% Perform MEX File check (also checks pre-reqs)
if (~mexFileCheck(localVer, cpath))
    return;
end

%Uninstall previous versions of OPTI
fprintf('\n- Checking for previous versions of OPTI Toolbox...\n');
no = opti_Uninstall('opti_Install.m',0);
if(no < 1)
    fprintf('Could not find a previous installation of OPTI Toolbox\n');
else
    fprintf('Successfully uninstalled previous version(s) of OPTI Toolbox\n');
end

%Add toolbox path to MATLAB
fprintf('\n- Adding OPTI Paths to MATLAB Search Path...');
genp = genpath(cd);
genp = regexp(genp,';','split');
%Folders to exclude from adding to Matlab path
i = 1;
rInd{:,:,i} = strfind(genp,'Documentation + Licenses'); i = i + 1;
rInd{:,:,i} = strfind(genp,'Source'); i = i + 1;
rInd{:,:,i} = strfind(genp,'CppAD'); i = i + 1;
rInd{:,:,i} = strfind(genp,'misdp'); i = i + 1;
rInd{:,:,i} = strfind(genp,'tex'); i = i + 1;
rInd{:,:,i} = strfind(genp,'.git'); i = i + 1;
rInd{:,:,i} = strfind(genp,'Crash Files'); i = i + 1;
rInd{:,:,i} = strfind(genp,'mexopts'); i = i + 1;
if(~exist([cd '\Solvers\Source\lib\win32\libclp.lib'],'file'))
    rInd{:,:,i} = strfind(genp,'Development'); i = i + 1;
end

ind = NaN(length(rInd{1}),1);
%Track indices of paths to remove from list
for i = 1:length(rInd{1})
    for j = 1:size(rInd,3)
        if(any(rInd{j}{i}))
            ind(i) = 1;
        end
    end
end

%Remove paths from above and add to matlab path
genp(ind == 1) = [];
addpath(genp{:});
rehash
fprintf('Done\n\n');
in = input('- Would You Like To Save the Path Changes? (Recommended) (y/n): ','s');
if(strcmpi(in,'y'))
    try
        savepath;
    catch %#ok<CTCH>
        warning('opti:install',['It appears you do not have administrator rights on your computer to save the Matlab path. '...
                                'In order to run OPTI Toolbox you will need to install it each time you wish to use it. To fix '...
                                'this please contact your system administrator to obtain administrator rights.']);
    end
end

%Post Install Test if requested
in = input('\n- Would You Like To Run Post Installation Tests? (Recommended) (y/n): ','s');
if(strcmpi(in,'y'))
    opti_Install_Test(1);
end

%Launch Examples page
web('https://inverseproblem.co.nz/OPTI/index.php/Examples/Examples','-browser');

%Finished
fprintf('\n\nOPTI Toolbox Installation Complete!\n');
disp('------------------------------------------------')

fprintf('\n\nYou now have the following solvers available to use:\n');
optiSolver;


function no = opti_Uninstall(token,del)
no = 0;
%Check nargin in, default don't delete and opti mode
if(nargin < 2 || isempty(del))
    del = 0;
end

%Check if we have anything to remove
paths = which(token,'-all');
len = length(paths);
%If mode is opti, should always be at least 1 if we are in correct directory
if(~len)
    error('Expected to find "%s" in the current directory - please ensure you are in the OPTI Toolbox directory');        
%If mode is opti, and there is one entry    
elseif(len == 1)
    %if len == 1, either we are in the correct folder with nothing to remove, or we are in the
    %wrong folder and there are files to remove, check CD
    if(any(strfind(paths{1},cd)))
        no = 0;
        return;
    else
        error('Expected to find "%s" in the current directory - please ensure you are in the OPTI Toolbox directory');
    end    
else %old ones to remove
    %Remove each folder found, and all subdirs under
    for n = 2:len
        %Absolute path to remove
        removeP = paths{n};
        %Search backwards for first file separator (we don't want the filename)
        for j = length(removeP):-1:1
            if(removeP(j) == filesep)
                break;
            end
        end
        removeP = removeP(1:max(j-1,1));        

        %Everything is lowercase to aid matching
        lrpath = lower(removeP);
        opath = regexp(lower(path),';','split');

        %Find & Remove Matching Paths
        no = 0;
        for i = 1:length(opath)
            %If we find it in the current path string, remove it
            fnd = strfind(opath{i},lrpath);        
            if(~isempty(fnd))  
                rmpath(opath{i});
                no = no + 1;
            end
        end
        
        %Check we aren't removing our development version
        rehash;
        if(isdir([removeP filesep 'Testing'])) %is this robust enough?
            fprintf('Found development version in "%s", skipping.\n',removeP);
            return;
        end

        %If delete is specified, also delete the directory
        if(del)
            stat = recycle; recycle('on'); %turn on recycling
            rmdir(removeP,'s'); %not sure if we dont have permissions here
            recycle(stat); %restore to original
        end
    end    
end

function matlabVerCheck()
mver = ver('MATLAB');
fprintf('\n- Checking MATLAB version and operating system...\n');
vv = regexp(mver.Version,'\.','split');
if(str2double(vv{1}) < 8)
    if(str2double(vv{2}) < 12)
        if(str2double(vv{2}) < 10)
            error('MATLAB 2011a or above is required to run OPTI - sorry!');
        else %2010a/b/sp1
            fprintf(2,'OPTI is designed for MATLAB 2011a or above.\nIt will install into 2010a, but you may experience reliability problems.\nPlease upgrade to R2011a or later.\n');
        end
    end
end

switch(mexext)
    case 'mexw32'
        error(['From v2.20 OPTI Toolbox only supports 64bit (Windows x64) platforms. Realistically, you should consider upgrading to 64bit for the best performance with OPTI.\n'...
               '\nIf however you would like to persist with 32bit, please download the last 32bit maintained version (v2.16) from the OPTI dropbox account:\n'...
               '%s'],'https://www.dropbox.com/s/ct2wmn1ajvujb3g/OptiToolbox_v2.16.zip?dl=0');
    case 'mexw64'
        fprintf('MATLAB %s 64bit (Windows x64) detected\n',mver.Release);
    otherwise
        error('OPTI Toolbox is compiled only for Windows systems - sorry!');
end


function OK = preReqChecks(cpath)
%Search for each required prereq
% Note we no longer search the registry, simply check if we can load a mex
% file which requires each runtime

if(~isempty(strfind(computer,'64')))
    arch = 'x64';
    icarch = 'Intel 64';
else
    arch = 'x86';
    icarch = 'IA32';
end

fprintf('\n- Checking for the required pre-requisites...\n');
missing = false;
cd('Utilities/Install');
[havVC,havIC,havIF] = opti_PreReqCheck(cpath);
cd(cpath);
%See if missing anything
if(~havVC || ~havIC || ~havIF)
    missing = true;
end

%Print Missing PreReqs
if(~havVC)
    fprintf(2,'Cannot find the Microsoft VC++ 2017 %s Redistributable!\n',arch); 
else
    fprintf('Found the Microsoft VC++ 2017 %s Redistributable\n',arch); 
end
% if(~havIC) %[not req from OPTI v >= 2.12]
%     fprintf(2,'Cannot find the Intel C++ XE 2013 %s Redistributable!\n',arch);
% else
%     fprintf('Found the Intel C++ XE 2013 %s Redistributable\n',arch); 
% end
if(~havIF)
    fprintf(2,'Cannot find the Intel Fortran XE 2017 %s Redistributable!\n',arch);
else
    fprintf('Found the Intel Fortran XE 2017 %s Redistributable\n',arch); 
end    

%Install Instructions for each Package
if(missing)
    fprintf(2,'\nYou are missing one or more pre-requisites. Please read the instructions below carefully to install them:\n\n');
    
    if(~havVC)
        fprintf(2,' Microsoft VC++ 2017:\n');
        switch(arch)
            case 'x64'
                fprintf(2,'- Download from: https://go.microsoft.com/fwlink/?LinkId=746572\n');
            case 'x86'
                fprintf(2,'- Download from: https://go.microsoft.com/fwlink/?LinkId=746571\n');
        end
        fprintf(2,['NOTE: If you have already downloaded and installed VC++ 2017 (and restarted MATLAB) - it may be that you are missing the Universal C Runtime (Universal CRT).\nThis is automatically installed '...
                    'with Windows Updates - but if you don''t have those turned on, you can download it from here:\nhttps://www.microsoft.com/en-us/download/details.aspx?id=48234\n\n']);
    end
    
%     if(~havIC) %[not req from OPTI v >= 2.12]
%         fprintf(2,' Intel C++ XE 2013:\n  - Download from: http://software.intel.com/en-us/articles/redistributable-libraries-for-intel-c-and-visual-fortran-composer-xe-2013-sp1-for-windows\n');
%         fprintf(2,'  - The download page will contain multiple links. Download the latest (highest number) update from the ''Intel C++ Composer XE 2013 for Windows Table''\n');
%         fprintf(2,'  - The download package will contain two files. Install the ''%s'' package.\n\n',icarch);
%     end
    
    if(~havIF) 
        fprintf(2,' Intel Fortran XE 2017:\n  - Download from: https://software.intel.com/en-us/articles/redistributables-for-intel-parallel-studio-xe-2017-composer-edition-for-windows\n');
        fprintf(2,'  - The download page will contain multiple links. Download the latest (highest number) update from the ''Intel Fortran Compiler for Windows Table''\n');
        fprintf(2,'  - The download package will contain two files. Install the ''%s'' package.\n\n',icarch);
    end
    
    fprintf(2,'\nOnce you have downloaded AND installed all the above packages, you MUST restart MATLAB.\n\nIf this message appears again after installing the above packages, try restarting your computer.\n\n\n');
    
    OK = false;
else
    OK = true;
end


function OK = mexFileCheck(localVer,cpath)

% Add paths required for checks
addpath([cd '/Solvers'])
addpath([cd '/Utilities'])
addpath([cd '/Utilities/opti'])

try
    % Check if we have the mex files (user may have cloned the repo without them)
    if (~exist(['clp.' mexext],'file') || ~exist(['ipopt.' mexext],'file') || ~exist(['asl.' mexext],'file'))
        fprintf('\n- New OPTI Installation Detected.\n');

        % We need to download the mex files / get the user to download them
        OK = downloadMexFiles(localVer);
        if (OK == false)
            return;
        end
    end
    
    % Do a pre-req check before attempting to read solver build versions
    OK = preReqChecks(cpath);
    if (OK == false)
        return; % can't check mex files if missing a pre req
    end

    fprintf('\n- Checking MEX File Release Information...\n');
    % Check if the current OPTI version matches the mex files
    [mexFilesOK, mexBuildVer] = checkMexFileVersion(localVer, false);
    if (mexFilesOK == false)
        % One or more mex files are out of date, report to the user
        if (isnan(mexBuildVer))
            fprintf(2,'One or more MEX files are not compatible with this version of OPTI\n');
        else
            fprintf(2,'One or more MEX files are not compatible with this version of OPTI (MEX v%.2f vs OPTI v%.2f)\n', mexBuildVer, localVer);
        end
        % See if the user wants to download them
        OK = downloadMexFiles(localVer);
    else
        % All up to date, nothing to check
        fprintf('MEX Files match OPTI Release.\n');
        OK = true;
    end

catch ME
    rmpath([cd '/Solvers']);
    rmpath([cd '/Utilities']);
    addpath([cd '/Utilities/opti'])
    rethrow(ME);
end


function [OK,buildVer] = checkMexFileVersion(localVer, verbose)

OK = true;
buildVer = localVer;
% Get a list of all solvers 
mexFiles = optiSolver('all');
% Add in utilities
mexFiles = [mexFiles,'asl','rmathlib'];

% Now search through and compare version info
for i = 1:length(mexFiles)
    if (~any(strcmpi(mexFiles{i},{'baron','cplex','matlab','gmatlab','mosek','sedumi'})))         
        % If SCIP the user may not have it
        if (strcmpi(mexFiles{i}, 'scip') && ~exist(['scip.' mexext], 'file'))
            continue;
        end
        try
            [~,optiBuildVer] = feval(lower(mexFiles{i}));
        catch 
            % Initially lots of mex files in the wild which don't support
            % the second arg, assume out of date
            OK = false;
            buildVer = NaN;
            if (verbose)
                fprintf('Error evaluating MEX File: %s\n', lower(mexFiles{i}));
            end
            return;
        end            
        if (optiBuildVer ~= localVer)
            if (optiBuildVer > localVer) % unusual case where user has newer mex files than opti source
                fprintf(2,'The MEX files you have appear to be for newer version of OPTI (MEX v%.2f vs OPTI v%.2f)\n', optiBuildVer, localVer);
                tellUserToUpdateOPTI();
            else    
                OK = false;
                buildVer = optiBuildVer;
                if (verbose)
                    fprintf('MEX File ''%s.%s'' is out of date (MEX v%.2f vs OPTI v%.2f)\n',mexFiles{i},mexext, optiBuildVer, localVer);
                else
                    break; % no point continuing check if not displaying
                end
            end
        end
    end
end


function OK = downloadMexFiles(localVer)

OK = true;
gitData = [];
% See if we can download directly from GitHub (2014b +)
fprintf('\n- Checking for updated MEX files from GitHub...');
if (exist('webread.m','file'))      
    try
        gitData = webread('https://api.github.com/repos/jonathancurrie/OPTI/releases/latest');
    catch
    end
end
if (isempty(gitData))
    % Cannot access internet / ML version too old, other error
    error('not implemented');
end

% If the Git version > local version, user needs to update OPTI source
if (gitVer > localVer)
    OK = false;
    tellUserToUpdateOPTI();
    return;
else  % download the latest files, even if local Ver > git Ver
    zipNameNoVer = ['optiMEXFiles_' mexext];
    numAssets = length(gitData.assets);
    mexFilesFoundOnGit = false;
    for i = 1:numAssets
        asset = gitData.assets(i);
        if (~isempty(asset))
            if (~isempty(strfind(asset.name, zipNameNoVer)))
                % Extract ver number from file name
                [~,fileName] = fileparts(asset.name);
                parts = regexp(fileName,'_','split');
                if (length(parts) == 4)
                    gitVer = str2double(parts{3}) + str2double(parts{4})/100;                                    
                    fprintf(' Found v%.2f\n', gitVer);
                else
                    % Should not happen...
                    fprintf(' Found\n');
                end
                % Start downloading
                fprintf('Downloading ''%s'' (%.2f MB), please wait...',asset.name,asset.size/(1024 * 1e3));
                tempLoc = [tempdir asset.name];
                try
                    websave(tempLoc,asset.browser_download_url);
                    mexFilesFoundOnGit = true;
                    fprintf(' Done!\n');
                catch ME
                    fprintf(' FAILED!\n');
                    fprintf(2, 'Failure: %s\n', ME.message);
                    fprintf(2, '\n\nPlease try running the installer again\n');
                    OK = false;
                    return;
                end                                   
                break;
            end
        end
    end
end

if (mexFilesFoundOnGit == true)
    % After a successful download, extract them
    [fdir,fname] = fileparts(tempLoc);
    unzipDir = [fdir filesep fname];
    % See if unzip directory already exists
    if (exist(unzipDir, 'dir'))
        % Delete folder and contents
        deleteFolder(unzipDir);
        rehash;
        rehash;
    end
    
    fprintf('Unzipping MEX Files...');
    unzip(tempLoc, unzipDir);
    rehash;
    rehash;
    fprintf(' Done!\n');
    
    % Now copy mex files
    try
        files = dir(unzipDir);
        for i = 1:length(files)
            file = files(i);        
            if (~file.isdir)
                fullFilePath = [file.folder filesep file.name];
                [~,mexName] = fileparts(file.name);
                copyTries = 0;
                ME = [];
                while (copyTries < 10) % annoyingly can fail now and then
                    copyTries = copyTries + 1;
                    try
                        clear(mexName);
                        if (any(strcmp(mexName,{'asl','coinR','coinW','mklJac','rmathlib'})))
                            destLoc = [cd filesep 'Utilities' filesep];
                        else
                            destLoc = [cd filesep 'Solvers' filesep];
                        end
                        if (movefile(fullFilePath, [destLoc file.name], 'f'))
                            break;
                        end
                    catch ME
                        rehash;
                        pause(0.001);
                    end
                end

                if (copyTries >= 10)
                    fprintf(2, 'Failed to copy MEX file ''%s'' from ''%s'' to ''%s''!\n\nEnsure you have adminstrator rights in your OPTI instllation directory.\n',...
                                file.name, fullFilePath, [destLoc file.name]);
                    if (~isempty(ME))
                        fprintf(2,'\n\nError: %s\n', ME.message);
                    end
                    OK = false;
                    break;
                end
            end
        end
        rehash;
        pause(0.01);
        rehash;        

        % Now delete unzip directory
        deleteFolder(unzipDir);        
        
    catch ME
        fprintf(2,'There was an error copying the MEX files. Please check the below error to see why it failed and correct as required:\n\n');
        fprintf(2,'Error: %s\n', ME.message);
        OK = false;
    end
else
    fprintf(2,'The OPTI MEX File package was not found in the latest release - please contact support. Sorry!\n');
    OK = false;
end


function OK = deleteFolder(delDir)

% Delete each file first
files = dir(delDir);
for i = 1:length(files)
    file = files(i);        
    if (~file.isdir)
        delete([file.folder filesep file.name]);        
    end
end
rehash;
rehash;

% Now delete the directory
OK = rmdir(delDir);
rehash;


