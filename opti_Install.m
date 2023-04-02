function opti_Install(savePath,runTests,openBrowser)
% OPTI Toolbox Installation File
%
%  opti_Install(savePath, runTests, openBrowser)
%
%   savePath: Save the paths added by OPTI to the MATLAB path 
%   runTests: Run the post-installation tests 
%   openBrowser: Whether to open the OPTI Toolbox Website after installation 
%
% All arguments are optional and if not supplied, the user will be prompted
% to enter their selection in the MATLAB Command Window. True is the
% default option for each argument.
%
% You MUST be in the current directory of this file!
%
%   Copyright (C) 2023 Jonathan Currie (Control Engineering)
%   https://controlengineering.co.nz/Wikis/OPTI/

% Handle missing input args
if (nargin < 3), openBrowser = []; end
if (nargin < 2), runTests = []; end
if (nargin < 1), savePath = []; end

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

%Uninstall previous versions of OPTI
fprintf('\n- Checking for previous versions of OPTI Toolbox...\n');
no = opti_Uninstall('opti_Install.m',0);
if(no < 1)
    fprintf('Could not find a previous installation of OPTI Toolbox\n');
else
    fprintf('Successfully uninstalled previous version(s) of OPTI Toolbox\n');
end

% Perform MEX File check (also checks pre-reqs)
if (~mexFileCheck(localVer, cpath))
    return;
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
if (isempty(savePath))
    in = input('- Would You Like To Save the Path Changes? (Recommended) (y/n): ','s');
else
    in = bool2yn(savePath);
end
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
if (isempty(runTests))
    in = input('\n- Would You Like To Run Post Installation Tests? (Recommended) (y/n): ','s');
else
    in = bool2yn(runTests);
end
if(strcmpi(in,'y'))
    opti_Install_Test(1);
end

%Launch Examples page
if (isempty(openBrowser) || (openBrowser == true))
    web('https://controlengineering.co.nz/Wikis/OPTI/index.php/Examples/Examples','-browser');
end

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

fprintf('\n- Checking MATLAB version and operating system...\n');
mver = ver('MATLAB');
% Sometimes we get multiple products here (no idea why), ensure we have MATLAB
if (length(mver) > 1)
    for i = 1:length(mver)
        if (strcmp(mver(i).Name, 'MATLAB'))
            mver = mver(i);
            break;
        end
    end
end

vv = regexp(mver.Version,'\.','split');
if(str2double(vv{1}) < 9) % https://au.mathworks.com/support/requirements/previous-releases.html
    fprintf(2,'OPTI is designed for MATLAB 2020b or above.\nIt will install into lower versions, but you may experience reliability problems.\nPlease upgrade to R2020b or later.\n');
end

switch(mexext)
    case 'mexw32'
        error('OPTI Toolbox is compiled only for 64bit systems - sorry!');
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
    fprintf(2,'Cannot find the Microsoft VC++ 2019 %s Redistributable!\n',arch); 
else
    fprintf('Found the Microsoft VC++ 2019 %s Redistributable\n',arch); 
end
% if(~havIC) %[not req from OPTI v >= 2.12]
%     fprintf(2,'Cannot find the Intel C++ XE 2013 %s Redistributable!\n',arch);
% else
%     fprintf('Found the Intel C++ XE 2013 %s Redistributable\n',arch); 
% end
if(~havIF)
    fprintf(2,'Cannot find the Intel Fortran XE 2019 %s Redistributable!\n',arch);
else
    fprintf('Found the Intel Fortran XE 2019 %s Redistributable\n',arch); 
end    

%Install Instructions for each Package
if(missing)
    fprintf(2,'\nYou are missing one or more pre-requisites. Please read the instructions below carefully to install them:\n\n');
    
    if(~havVC)
        fprintf(2,' Microsoft VC++ 2019:\n');
        switch(arch)
            case 'x64'
                fprintf(2,'- Download from: https://aka.ms/vs/17/release/vc_redist.x64.exe\n');
            case 'x86'
                fprintf(2,'- Download from: https://aka.ms/vs/17/release/vc_redist.x86.exe\n');
        end
        fprintf(2,['NOTE: If you have already downloaded and installed VC++ 2019 (and restarted MATLAB) - it may be that you are missing the Universal C Runtime (Universal CRT).\nThis is automatically installed '...
                    'with Windows Updates - but if you don''t have those turned on, you can download it from here:\nhttps://www.microsoft.com/en-us/download/details.aspx?id=48234\n\n']);
    end
    
%     if(~havIC) %[not req from OPTI v >= 2.12]
%         fprintf(2,' Intel C++ XE 2013:\n  - Download from: http://software.intel.com/en-us/articles/redistributable-libraries-for-intel-c-and-visual-fortran-composer-xe-2013-sp1-for-windows\n');
%         fprintf(2,'  - The download page will contain multiple links. Download the latest (highest number) update from the ''Intel C++ Composer XE 2013 for Windows Table''\n');
%         fprintf(2,'  - The download package will contain two files. Install the ''%s'' package.\n\n',icarch);
%     end
    
    if(~havIF) 
        fprintf(2,' Intel Fortran XE 2019:\n  - Download from: https://software.intel.com/en-us/articles/redistributable-libraries-for-intel-c-and-fortran-2019-compilers-for-windows\n');
        fprintf(2,'  - The download page will contain multiple links. Download the latest (highest number) update from the ''Intel Fortran Compiler for Windows Table''\n');
        fprintf(2,'  - The download package will contain two files. Install the ''%s'' package.\n\n',icarch);
    end
    
    fprintf(2,'\nOnce you have downloaded AND installed all the above packages, you MUST restart MATLAB.\n\nIf this message appears again after installing the above packages, try restarting your computer.\n\n\n');
    
    OK = false;
else
    OK = true;
end


function OK = mexFileCheck(localVer,cpath)

% Check if a dev version of OPTI
if (exist([cd '/Solvers/Source/lib/win64/libclp.lib'],'file'))
    fprintf('\nOPTI Development Version Detected, Skipping MEX File Check\n');
    OK = true;
    return;
end

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
                fprintf(2,'The MEX files you have are for newer version of OPTI (MEX v%.2f vs Local OPTI v%.2f)\n', optiBuildVer, localVer);
                fprintf(2, 'You will need to update your version of OPTI, please follow the instructions below:\n');
                opti_printUpdateInfo();
                error('OPTI Install Error: OPTI source must be updated');
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

function r = optiRound(r, n)
mver = ver('MATLAB');
vv = regexp(mver.Version,'\.','split');
needOptiRound = false;
% MATLAB 2014b  introduced round(r,n), before then approximate it
if(str2double(vv{1}) < 8)
    needOptiRound = true;
elseif (str2double(vv{1}) == 8 && str2double(vv{2}) < 4)
    needOptiRound = true;
end
if (needOptiRound == true)
    r = round(r*10^n)/(10^n);
else
    r = round(r,n);
end


function OK = downloadMexFiles(localVer)

localVer = optiRound(localVer, 3);
gitData = [];
mexFilesFoundOnGit = false;
% See if we can download directly from GitHub (2014b +)
if (exist('webread.m','file'))  
    % See if the user wants us to automatically download the mex files
    in = input('\n- Would You Like OPTI To Attempt to Download the MEX Files Automatically? (Recommended) (y/n): ','s');
    if (strcmpi(in,'y'))
        fprintf('\n- Checking for updated MEX files from GitHub...');
        try
            gitData = webread('https://api.github.com/repos/jonathancurrie/OPTI/releases/latest');
        catch
        end
    end
end

% If we got the download info from Github
if (~isempty(gitData))
    % Download the latest files
    zipNameNoVer = ['optiMEXFiles_' mexext];
    numAssets = length(gitData.assets);
    for i = 1:numAssets
        asset = gitData.assets(i);
        if (~isempty(asset))
            if (~isempty(strfind(asset.name, zipNameNoVer)))
                % Extract ver number from file name
                [~,fileName] = fileparts(asset.name);
                parts = regexp(fileName,'_','split');
                if (length(parts) == 4)
                    gitVer = optiRound(str2double(parts{3}) + str2double(parts{4})/100, 3);                                    
                    fprintf(' Found v%.2f\n', gitVer);
                    % If the Git version > local version, user needs to update OPTI source
                    if (gitVer > localVer)
                        OK = false;
                        fprintf(2, 'The GitHub MEX Files are for a newer version of OPTI (GitHub %.2f, Local OPTI %.2f)\n', gitVer, localVer);
                        fprintf(2, 'You will need to update your version of OPTI, please follow the instructions below:\n');
                        opti_printUpdateInfo();
                        return;
                    elseif(gitVer < localVer) % should not happen
                        fprintf(2, 'Your version of OPTI (%.2f) is newer than the version of MEX files available (%.2f).\n', localVer, gitVer);
                        fprintf(2, 'This normally occurs if you are working on the develop branch.\n');
                        fprintf(2, 'OPTI will download the latest release version, but please check GitHub regularly for an updated release.\n\n');
                    end
                else
                    % Should not happen...
                    fprintf(' Found\n');
                end
                % Start downloading
                fprintf('Downloading ''%s'' (%.2f MB), this will take a few minutes, please wait...',asset.name,asset.size/(1024 * 1e3));
                tempLoc = [tempdir asset.name];
                try
                    websave(tempLoc,asset.browser_download_url);
                    mexFilesFoundOnGit = true;
                    fprintf(' Done!\n');
                catch ME
                    fprintf(' FAILED!\n');
                    fprintf(2, 'Failure: %s\n', ME.message);
                    fprintf(2, '\n\nPlease follow the below instructions to download the MEX files manually.\n');
                end                                   
                break;
            end
        end
    end
end

if (mexFilesFoundOnGit == true)
    OK = copyMexFiles(tempLoc);
else
    % Cannot access internet / ML version too old, other error
    fprintf('\nIn order to update the MEX files in your OPTI installation, please visit:\n');
    disp(' <a href="https://github.com/jonathancurrie/OPTI/releases/latest">https://github.com/jonathancurrie/OPTI/releases/latest</a>')
    fprintf('and download "optiMEXFiles_%s_x_xx.zip" where x_xx is the latest release number.\n', mexext);
    input('\nOnce you have downloaded the zipped files, press enter to show OPTI where they are located (press enter to continue):  ', 's'); 
    [FileName,PathName] = uigetfile('*.zip', 'Select the OPTI MEX Files Download');

    % Check seems a reasonable folder for opti files
    if (~ischar(FileName))
        error('OPTI cannot continue without knowing where you have downloaded the MEX files to. Please re-run the installer to continue.');
    end
    [~,~,ext] = fileparts(FileName);
    if (isempty(strfind(FileName,['optiMEXFiles_' mexext])) || ~strcmp(ext, '.zip'))
        error('OPTI did not recognise ''%s'' as a valid OPTI MEX File Download. Please run the installer again to select the correct download.',FileName);
    end
    % Copy them over
    OK = copyMexFiles([PathName FileName]);
end


function OK = copyMexFiles(loc)

OK = true;

% After a successful download, extract them
[fdir,fname] = fileparts(loc);
unzipDir = [fdir filesep fname];
% See if unzip directory already exists
if (exist(unzipDir, 'dir'))
    % Delete folder and contents
    deleteFolder(unzipDir);
    rehash;
    rehash;
end

fprintf('Unzipping MEX Files (please wait)...');
unzip(loc, unzipDir);
rehash;
rehash;
fprintf(' Done!\n');

% Now copy mex files
try
    files = dir(unzipDir);
    for i = 1:length(files)
        file = files(i);        
        if (~file.isdir)
            if (~isfield(file,'folder'))
                fullFilePath = [unzipDir filesep file.name];
            else
                fullFilePath = [file.folder filesep file.name];
            end
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
    
    if (OK)
        rehash;
        pause(0.01);
        rehash;        

        % Now delete unzip directory
        OK = deleteFolder(unzipDir);        
    end

catch ME
    fprintf(2,'There was an error copying the MEX files. Please check the below error to see why it failed and correct as required:\n\n');
    fprintf(2,'Error: %s\n', ME.message);
    OK = false;
end


function OK = deleteFolder(delDir)

% Delete each file first
files = dir(delDir);
for i = 1:length(files)
    file = files(i);        
    if (~file.isdir)
        if (~isfield(file,'folder'))
            delete([delDir filesep file.name]);        
        else
            delete([file.folder filesep file.name]);        
        end
    end
end
rehash;
rehash;

% Now delete the directory
OK = rmdir(delDir);
rehash;


function in = bool2yn(val)
if (isempty(val) || val == true)
    in = 'y';
else
    in = 'n';
end
