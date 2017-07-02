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
cur_ver = optiver();

fprintf('\n------------------------------------------------\n')
fprintf(['  INSTALLING OPTI TOOLBOX ver ' sprintf('%1.2f',cur_ver) '\n'])

cd(cpath);
% Check if we need to download the MEX files
if (~exist([cd filesep 'Solvers/Source/lib/win64/libclp.lib'],'file')) % skip on dev machine
    fprintf('\n- Checking for updated MEX files from GitHub...\n');
    try
        if (exist('webread.m','file'))      
            try
                gitData = webread('https://api.github.com/repos/jonathancurrie/OPTI/releases/latest');
            catch ME
                fprintf(2,'There was an error querying GitHub for the latest MEX file information. Please ensure you are connected to the internet!\n');
                rethrow(ME);
            end
        else
            fprintf(2,'Your version of MATLAB does not support the required OPTI GitHub interface. Please update your MATLAB version.\n');
            error('OPTI Install Error: MATLAB Version does not support webread()');
        end
        if (isempty(gitData))
            error('OPTI Install Error: No valid data received from GitHub!');
        end
        nameComp = regexp(gitData.name,' ','split');
        if (length(nameComp) ~= 3)
            error('OPTI Install Error: The latest version name (%s) is not compatible! Please report this error.',gitData.name);
        end
        verNum = str2double(nameComp{3});
        if (isnan(verNum))
            verNum = str2double(nameComp{3}(2:end));
            if (isnan(verNum))
                error('OPTI Install Error: Could not convert latest release version number (%s) to double! Please report this error.',nameComp{3});
            end
        end
        doDownload = false;
        if (~exist(['rmathlib.' mexext], 'file'))
            if (verNum ~= cur_ver)
                fprintf(2,'Your version of OPTI is not the most recent, please update it from GitHub (https://github.com/jonathancurrie/OPTI) before continuing.\n');
                error('OPTI Install Error: OPTI version mismatch');
            else
                fprintf('A new OPTI installation is detected, the required solver MEX files will be automatically downloaded from GitHub:\n');
                doDownload = true;
            end
        elseif(checkRMathlibBuildDate(gitData)) % compare the build date in the mex file
            doDownload = true;
        end
        % Download as required
        if (doDownload)
            downloadMEXFiles(gitData);
        end
    catch ME
        fprintf(2,'If you continue to see this error, you may manually download the required MEX files (all) from:\nhttps://github.com/jonathancurrie/OPTI/releases/latest\n\n\n');
        rethrow(ME);
    end
end
%Perform pre-req check
if(~preReqChecks(cpath))
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

fprintf('\n- Checking for required pre-requisites...\n');
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



function needUpdate = checkRMathlibBuildDate(gitData)

needUpdate = false;
foundAsset = false;
% Find ASL asset
numAssets = length(gitData.assets);
for i = 1:numAssets   
    asset = gitData.assets(i);
    if(~isempty(asset))
        [~,name,ext] = fileparts(asset.name);
        if (strcmp(name, 'rmathlib'))
            foundAsset = true;
            clear('rmathlibTemp');
            websave([cd filesep 'Utilities' filesep 'rmathlibTemp' ext],asset.browser_download_url); % grab it
            rehash;
            rehash;
            try
                [~,localTS] = rmathlib;
                [~,refTS] = rmathlibTemp;                
                refDateTime = datenum(refTS);
                localDateTime = datenum(localTS);
                if (refDateTime > localDateTime)
                    needUpdate = true;
                end
            catch
                needUpdate = true; % safe option - not everyone will have the updated mex file for a while
            end
        end
    end
end
if (foundAsset == false)
    % Failed by here
    error('OPTI Install Error: Could not find rmathlib in the latest release information! Please report this error');
else
    try
        delete(which(['rmathlibTemp.' mexext])); 
    catch
    end
end


function downloadMEXFiles(gitData)

% For each MEX file
numAssets = length(gitData.assets);
for i = 1:numAssets   
    asset = gitData.assets(i);
    if(~isempty(asset))
        fprintf('Downloading %d of %d %18s (%.2f MB)...',i,numAssets,asset.name,asset.size/(1024 * 1e3));
        [~,name] = fileparts(asset.name);
        % See if utility or solver
        isUtil = false;
        if (any(strcmp(name,{'asl','coinR','coinW','mklJac','rmathlib'})))
            isUtil = true;
        end
        clear(name);
        if (isUtil)
            websave([cd filesep 'Utilities' filesep asset.name],asset.browser_download_url);
        else
            websave([cd filesep 'Solvers' filesep asset.name],asset.browser_download_url);
        end
        fprintf(' Done!\n');
    end
end
    
