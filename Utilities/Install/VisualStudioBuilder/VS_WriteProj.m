function [projPath,guid] = VS_WriteProj(srcpath,projName,incpath,opts)
%WRITEVS  Create a Visual Studio Project from selected paths
%
%   This function attempts to automatically create a Visual Studio Project
%   from a supplied source directory path (or paths). It works for me, but
%   may not for every project! Use with caution. Project will be created
%   one directory up from srcpath.
%
%   proj_path = VS_WriteProj(source_path,proj_name,include_path,options)
%
%   source_path     - absolute path to the source code to compile (may be a
%                     cell array of locations).
%   proj_name       - project name (optional).
%   include_path    - absolute path of directory to add to C/C++ Inlude
%                     Directories (may be a cell array of locations). Note
%                     files in these directories are not added to the
%                     project tree. (optional)
%   options         - structue with fields: (optional)
%                       'cpp'      read C++ and C files (true), read FORTRAN files (false) {true}
%                       'exPP'     cell array of extra preprocessor definitions {[]}
%                       'ex64PP'   cell array of extra preprocessor definitions for 64 bit builds (in addition to exPP) {[]}
%                       'openMP'   true / false {false}
%                       'mkl'      true / false to include MKL headers {false}
%                       'charset'  Character Set {'unicode'}, 'multibyte'
%                       'toolset'  'v100' for VS2010, 'v110' for VS2012, v120 for VS2013, v140 for VS2015, v141 for VS2017 {v141}
%                       'ifortver' Intel Fortran version 'XE11' for XE2011/SP1, 'XE13' for XE2013/SP1, 'XE15' for XE2015, 'XE16' for XE2016 {XE16}
%                       'exclude'  cell array of source files to exclude from project {[]}
%                       'exFilter' cell array of filtered source files (e.g. abc*) to exclude from project {[]}
%                       'exFolder' cell array of source folders to exclude from project {[]}
%                       'include'  cell array of source files to include in project (only these will be added) {[]}                     
%                       'console'  true / false {false}
%                       'linklib'  cell array of libraries to add to linker line {[]}
%                       'linkpath' cell array of directories to add to linker search path {[]}
%                       'mkllink'  true / false to link MKL {false}
%                       'outname'  Output file name (rather than proj_name above)
%                       'dll'      true / false to build as a DLL {false}
%                       'def'      module definition file (only compatible with DLL projects) {[]}
%                       'compileAsCpp' compile .c as cpp
%                       'ExcepwCExtern' Accept C exceptions as extern

%Extract first source path if multiple paths specified (first is assumed base)
allpaths = srcpath;
if(iscell(srcpath))
    if(size(allpaths,2) > 1)
        allpaths = allpaths';
    end
    srcpath = allpaths{1};
end

%Error checking + default args
if(isempty(srcpath) || ~exist(srcpath,'dir')), error('Cannot find the directory:\n%s',srcpath); end
if(nargin < 2 || isempty(projName)), projName = '.lib'; end
if(nargin < 3), incpath = []; end
if(~isfield(opts,'cpp')), opts.cpp = true; end
if(~isfield(opts,'exPP')), opts.exPP = []; end
if(~isfield(opts,'ex64PP')), opts.ex64PP = []; end
if(~isfield(opts,'openMP')), opts.openMP = false; end
if(~isfield(opts,'charset')), opts.charset = 'unicode'; end
if(~isfield(opts,'mkl')), opts.mkl = false; end
if(~isfield(opts,'toolset')), opts.toolset = 'v141'; end
if(~isfield(opts,'ifortver')), opts.ifortver = 'XE16'; end
if(~isfield(opts,'include')), opts.include = []; end
if(~isfield(opts,'exclude')), opts.exclude = []; end
if(~isfield(opts,'exFolder')), opts.exFolder = []; end
if(~isfield(opts,'exFilter')), opts.exFilter = []; end
if(~isfield(opts,'console')), opts.console = false; end
if(~isfield(opts,'linklib')), opts.linklib = []; end
if(~isfield(opts,'linkpath')), opts.linkpath = []; end
if(~isfield(opts,'mkllink')), opts.mkllink = false; end
if(~isfield(opts,'outname')), opts.outname = projName; end
if(~isfield(opts,'dll')), opts.dll = false; end
if(~isfield(opts,'def')), opts.def = []; end
if(~isfield(opts,'empty')), opts.empty = false; end
if(~isfield(opts,'archname')), opts.archname = false; end
if(~isfield(opts,'compileAsCpp')), opts.compileAsCpp = false; end
if(~isfield(opts,'ExcepwCExtern')), opts.ExcepwCExtern = false; end

%If VS2010 (vs100), remove from options (not required)
if(strcmpi(opts.toolset,'v100')), opts.toolset = []; end
%Check charset
switch(lower(opts.charset))
    case 'unicode'
        opts.charset = 'Unicode';
    case 'multibyte'
        opts.charset = 'MultiByte';
    otherwise
        error('Unknown character set selected, valid options are Unicode or Multibyte');
end
%Setup configuration type
if(opts.console)
    opts.configType = 'Application';
elseif(opts.dll)
    opts.configType = 'DynamicLibrary';
else
    opts.configType = 'StaticLibrary';
end
%Check cell arguments
if(~isempty(opts.linklib))
    if(~iscell(opts.linklib)), opts.linklib = {opts.linklib}; end
    if(size(opts.linklib,2) > 1), opts.linklib = opts.linklib'; end
end
if(~isempty(opts.linkpath))
    if(~iscell(opts.linkpath)), opts.linkpath = {opts.linkpath}; end
    if(size(opts.linkpath,2) > 1), opts.linkpath = opts.linkpath'; end
end
if(~isempty(opts.include))
    if(~iscell(opts.include)), opts.include = {opts.include}; end
    if(size(opts.include,2) > 1), opts.include = opts.include'; end
end

%Do the same for extra include paths (paths not added to the project tree)
incpaths = incpath;
if(iscell(incpath))
    if(size(incpaths,2) > 1)
        incpaths = incpaths';
    end
end

%Remove trailing \
if(iscell(allpaths))
    for i = 1:size(allpaths,1)
        if(allpaths{i}(end) == '\')
            allpaths{i}(end) = [];
        end
    end
else
    if(srcpath(end) == '\')
        srcpath(end) = [];
        allpaths(end) = [];
    end    
end
if(iscell(incpaths))
    for i = 1:size(incpaths,1)
        if(incpaths{i}(end) == '\')
            incpaths{i}(end) = [];
        end
    end
elseif(~isempty(incpaths))
    if(incpaths(end) == '\')
        incpaths(end) = [];
    end    
end

%Create Project Directory
if(any(strfind(srcpath,':'))) %assume full path
    %Go up one folder by removing last folder line
    ipath = srcpath;
    if(~opts.empty)
        ind = strfind(srcpath,'\');
        ipath(ind(end):end) = [];
    end
    %Make the directory
    mkdir(ipath,projName);
    projPath = [ipath '\' projName];
else    
    mkdir(cd,projName);    
    projPath = [cd '\' projName];
end

%IF empty project AND DLL, create a dummy main
if(opts.dll && opts.empty)
    try
        fp = fopen([ipath filesep projName filesep 'main.cpp'],'w');
        fprintf(fp,'int main() {\n  return 0;\n}\n');
        fclose(fp);
    catch
        fclose(fp);
    end
    allpaths = [ipath '\' projName];
end

% THE REMAINDER OF THIS FILE IMPLEMENTS A STATIC LIBRARY WITH 32BIT AND 64BIT CONFIGURATIONS. 
% IT MAY NOT WORK WITH OLDER (OR FUTURE) VERSIONS OF VISUAL STUDIO

%Get a Project GUID
guid = getProjGUID();

%Header
if(opts.cpp)
    docNode = com.mathworks.xml.XMLUtils.createDocument('Project');
    p = docNode.getDocumentElement;
    p.setAttribute('DefaultTargets','Build');
    if(strcmpi(opts.toolset,'v141'))
        p.setAttribute('ToolsVersion','15.0');
    elseif(strcmpi(opts.toolset,'v140'))
        p.setAttribute('ToolsVersion','14.0');
    elseif(strcmpi(opts.toolset,'v120'))
        p.setAttribute('ToolsVersion','12.0');
    else
        p.setAttribute('ToolsVersion','4.0');
    end
    p.setAttribute('xmlns','http://schemas.microsoft.com/developer/msbuild/2003');
    %Project Configuration
    pc = createSection(docNode,'ItemGroup','ProjectConfigurations');
    pc.appendChild(writeProjConfig(docNode,'Debug','Win32'));
    pc.appendChild(writeProjConfig(docNode,'Debug','x64'));
    pc.appendChild(writeProjConfig(docNode,'Release','Win32'));
    pc.appendChild(writeProjConfig(docNode,'Release','x64'));
    p.appendChild(pc);    
    %Globals
    pc = createSection(docNode,'PropertyGroup','Globals');
    addElemText(docNode,pc,'ProjectGuid',['{' guid '}']); 
    addElemText(docNode,pc,'Keyword','Win32Proj');
    addElemText(docNode,pc,'RootNamespace',projName);
    p.appendChild(pc);
    %Import
    p.appendChild(createImport(docNode,'$(VCTargetsPath)\Microsoft.Cpp.Default.props',[],[]));        
    %Config Setup
    debug = struct('debug',true,'optimize',false);
    release = struct('debug',false,'optimize',true);
    p.appendChild(createConfig(docNode,'Win32',opts,debug));
    p.appendChild(createConfig(docNode,'x64',opts,debug));
    p.appendChild(createConfig(docNode,'Win32',opts,release));
    p.appendChild(createConfig(docNode,'x64',opts,release));

    %Imports
    p.appendChild(createImport(docNode,'$(VCTargetsPath)\Microsoft.Cpp.props',[],[]));
    pc = docNode.createElement('ImportGroup');
    pc.setAttribute('Label','ExtensionSettings');
    p.appendChild(pc);
    p.appendChild(createPropSheet(docNode,'Debug','Win32'));
    p.appendChild(createPropSheet(docNode,'Debug','x64'));
    p.appendChild(createPropSheet(docNode,'Release','Win32'));
    p.appendChild(createPropSheet(docNode,'Release','x64'));

    %Macros
    pc = docNode.createElement('PropertyGroup');
    pc.setAttribute('Label','UserMacros');
    p.appendChild(pc);
    %Custom Name
    if(~isempty(opts.outname))
        configs = {'Debug|Win32','Debug|x64','Release|Win32','Release|x64'};
        for i = 1:length(configs)
            pc = createSection(docNode,'PropertyGroup','Condition',['''$(Configuration)|$(Platform)''==''' configs{i} '''']);
            if(opts.archname)
                addElemText(docNode,pc,'TargetName',[opts.outname configs{i}(end-1:end)]);  
            else
                addElemText(docNode,pc,'TargetName',opts.outname);        
            end
            p.appendChild(pc);
        end
    end
else %FORTRAN
    docNode = com.mathworks.xml.XMLUtils.createDocument('VisualStudioProject');
    p = docNode.getDocumentElement;
    p.setAttribute('ProjectType','typeStaticLibrary');
    p.setAttribute('ProjectCreator','Intel Fortran');
    if(strcmpi(opts.ifortver,'XE15'))
        p.setAttribute('Keyword','Static Library');
    else
        p.setAttribute('Keyword','StaticLibrary');
    end
    p.setAttribute('Version','11.0');
    p.setAttribute('ProjectIdGuid',['{' guid '}']);
    pc = createSection(docNode,'Platforms');
    pc.appendChild(createSection(docNode,'Platform','Name','Win32'));
    pc.appendChild(createSection(docNode,'Platform','Name','x64'));
    p.appendChild(pc);
end

%Read in Source + Header Files from Main Path(s)
% if(~opts.empty)
    [src,hdr] = VS_BuildFileList(allpaths,opts.exFolder,opts.cpp);
    if(~isempty(opts.include))
        %Only source files specified will be added
        src = {allpaths, opts.include, 1}; 
    else
        %Remove files in exclude list
        if(~isempty(opts.exclude))
            for i = 1:size(src,1)
                if(~isempty(src{i,2}))
                    ind = ismember(src{i,2},opts.exclude);
                    if(any(ind))
                        src{i,2}(ind) = [];
                    end
                end
            end
        end
        %Remove files in exclude filter list
        if(~isempty(opts.exFilter))
            for i = 1:size(src,1)
                for j = 1:length(opts.exFilter)
                    if(~isempty(src{i,2}))
                        ind = regexp(src{i,2},opts.exFilter{j});
                        lind = true(size(ind));
                        for k = 1:length(ind)
                            if(ind{k} == 1)
                                lind(k) = false;
                            end
                        end
                        src{i,2} = src{i,2}(lind);
                    end
                end
            end
        end
    end
    
    %Read in Extra Header Files to Add
    projhdr = hdr;
    [~,hdrinc] = VS_BuildFileList(incpaths,opts.exFolder,opts.cpp);
    hdr = [hdr;hdrinc];
% else
%     src = [];
%     hdr = [];
%     projhdr = [];
% end


%Add MKL include path if requested
if(opts.mkl || opts.mkllink)
    [~,~,mkl_inc,mkl_lib,mkl_cmplr] = opti_FindMKL();
    if(isempty(hdr))
        hdr = {mkl_inc {'mkl.h'} 1};
    else
        no = hdr{end,end}+1;
        hdr = [hdr;{mkl_inc {'mkl.h'} no}];
    end
    if(opts.mkllink)
        mkl32libs = {'mkl_intel_c.lib', 'mkl_intel_thread.lib', 'mkl_core.lib',  'libiomp5md.lib'}';
        mkl64libs = {'mkl_intel_lp64.lib', 'mkl_intel_thread.lib', 'mkl_core.lib',  'libiomp5md.lib'}';
        switch(computer)
            case 'PCWIN'
                mkl32paths = {mkl_lib,mkl_cmplr}';
                mkl64paths = {regexprep(mkl_lib,'ia32','intel64'),regexprep(mkl_cmplr,'ia32','intel64')}';
            case 'PCWIN64'
                mkl64paths = {mkl_lib,mkl_cmplr}';
                mkl32paths = {regexprep(mkl_lib,'intel64','ia32'),regexprep(mkl_cmplr,'intel64','ia32')}';
        end
        
        if(isempty(opts.linklib))
            opts.link32lib = mkl32libs;
            opts.link64lib = mkl64libs;
        else
            opts.link32lib = [opts.linklib; mkl32libs];
            opts.link64lib = [opts.linklib; mkl64libs];
        end
        if(isempty(opts.linkpath))
            opts.link32path = mkl32paths;
            opts.link64path = mkl64paths;
        else
            opts.link32path = [opts.linkpath; mkl32paths];
            opts.link64path = [opts.linkpath; mkl64paths];
        end
    else
        opts.link32lib = opts.linklib;
        opts.link64lib = opts.linklib;
        opts.link32path = opts.linkpath;
        opts.link64path = opts.linkpath;
    end
else
    opts.link32lib = opts.linklib;
    opts.link64lib = opts.linklib;
    opts.link32path = opts.linkpath;
    opts.link64path = opts.linkpath;
end

if(opts.archname && ~isempty(opts.def))
    opts.def32 = regexprep(opts.def,'.def','32.def');
    opts.def64 = regexprep(opts.def,'.def','64.def');
else
    opts.def32 = opts.def;
    opts.def64 = opts.def;
end

%Combine pp for win32 and win64
if(isempty(opts.exPP))
    opts.ex32PP = [];    
elseif(~isempty(opts.exPP) && isempty(opts.ex64PP))
    opts.ex32PP = opts.exPP;
    opts.ex64PP = opts.exPP;
elseif(~isempty(opts.exPP) && ~isempty(opts.ex64PP))
    opts.ex32PP = opts.exPP;
    if(~iscell(opts.exPP)), opts.exPP = {opts.exPP}; end
    if(~iscell(opts.ex64PP)), opts.ex64PP = {opts.ex64PP}; end
    opts.ex64PP = [opts.exPP opts.ex64PP];
end 

if(opts.cpp)
    %Debug Detailed Settings
    p.appendChild(writeDebugDetail(docNode,allpaths,hdr,'Win32',opts.ex32PP,opts.openMP,opts.console,opts.link32lib,opts.link32path,opts.def32,opts.compileAsCpp,opts.ExcepwCExtern));
    p.appendChild(writeDebugDetail(docNode,allpaths,hdr,'x64',opts.ex64PP,opts.openMP,opts.console,opts.link64lib,opts.link64path,opts.def64,opts.compileAsCpp,opts.ExcepwCExtern));
    %Release Detailed Settings
    p.appendChild(writeReleaseDetail(docNode,allpaths,hdr,'Win32',opts.ex32PP,opts.openMP,opts.console,opts.link32lib,opts.link32path,opts.def32,opts.compileAsCpp,opts.ExcepwCExtern));
    p.appendChild(writeReleaseDetail(docNode,allpaths,hdr,'x64',opts.ex64PP,opts.openMP,opts.console,opts.link64lib,opts.link64path,opts.def64,opts.compileAsCpp,opts.ExcepwCExtern));

    %Write VS Filters
    VS_WriteFilters(projPath,projName,allpaths,src,projhdr);

    %Source Files
    p.appendChild(createFileList(docNode,allpaths,src,'ClCompile'));
    %Header Files
    p.appendChild(createFileList(docNode,allpaths,hdr,'ClInclude'));

    %Imports
    p.appendChild(createImport(docNode,'$(VCTargetsPath)\Microsoft.Cpp.targets',[],[]));
    pc = docNode.createElement('ImportGroup');
    pc.setAttribute('Label','ExtensionTargets');
    p.appendChild(pc);
else %FORTRAN
    %Configurations
    pc = createSection(docNode,'Configurations');
    pc.appendChild(writeIFortConfig(docNode,'Debug','Win32',opts.ex32PP,incpaths,opts.ifortver));
    pc.appendChild(writeIFortConfig(docNode,'Release','Win32',opts.ex64PP,incpaths,opts.ifortver));
    pc.appendChild(writeIFortConfig(docNode,'Debug','x64',opts.ex32PP,incpaths,opts.ifortver));
    pc.appendChild(writeIFortConfig(docNode,'Release','x64',opts.ex64PP,incpaths,opts.ifortver));
    p.appendChild(pc);
    %Write VS Filters
%     VS_WriteFilters(projPath,projName,allpaths,src,projhdr);
    %Files
    pc = createSection(docNode,'Files');
    pc.appendChild(createMulSection(docNode,'Filter',{{'Name','Header Files'},{'Filter','fi;fd'}}));
    pc.appendChild(createMulSection(docNode,'Filter',{{'Name','Resource Files'},{'Filter','rc;ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe'}}));
    pc.appendChild(createMulSection(docNode,'Filter',{{'Name','Source Files'},{'Filter','f90;for;f;fpp;ftn;def;odl;idl;F'}}));
    createFileList(docNode,allpaths,src,'File',true,pc);
%     createFileList(docNode,allpaths,hdr,'File',true,pc);
    p.appendChild(pc);
end

%Write Whole Project
xmlwrite([projPath '\' projName '.xml'],docNode);
if(opts.cpp)
    xmlwrite([projPath '\' projName '.vcxproj'],docNode);
else
    xmlwrite([projPath '\' projName '.vfproj'],docNode);
end

%Create File List Section (lists files to include / compile)
function pc = createFileList(docNode,upath,files,elem,isFort,pRoot)
if(nargin < 5), isFort = false; end
if(nargin < 6), pRoot = false; end
if(~isFort)
    pc = docNode.createElement('ItemGroup');
end
for i = 1:size(files,1) %for each path
    if(iscell(upath))
        p = upath{files{i,3}};
    else
        p = upath;
    end

    %Check we are on the same path as our base one
    h = files{i,1};
    if(isempty(strfind(h,p))) %empty means not same base bath
        attp = getPathSep(h,p);
    else %OK contained with base path
        if(any(strfind(p,':')))
            ind = strfind(p,'\'); %find end file separator
            len2 = ind(end)+1; %drop last folder (we are always back one)
        else
            len2 = 1;
        end    
        attp = files{i,1}(len2:end);
    end
    if(~isempty(files{i,2}))
        for j = 1:size(files{i,2},1) %for each file
            if(~isFort)
                pc1 = docNode.createElement(elem);
                pc1.setAttribute('Include',['..\' attp '\' files{i,2}{j}]);
                pc.appendChild(pc1);
            else
                pc = docNode.createElement(elem);
                pc.setAttribute('RelativePath',['..\' attp '\' files{i,2}{j}]);
                pRoot.appendChild(pc);
            end
        end
    end
end

%Create an XML Section
function pc = createSection(docNode,name,labelname,labelval)
pc = docNode.createElement(name);
switch(nargin)
    case 3
        pc.setAttribute('Label',labelname);
    case 4
        pc.setAttribute(labelname,labelval);
end

%Create an XML Section with multiple attributes
function pc = createMulSection(docNode,name,attr)
pc = docNode.createElement(name);
for i = 1:length(attr)
    pc.setAttribute(attr{i}{1},attr{i}{2});
end

%Lowlevel routine to add an element with text
function addElemText(docNode,pc,elem,text)
c = docNode.createElement(elem);
c.appendChild(docNode.createTextNode(text));
pc.appendChild(c);

%Create Build Configuration Section
function pc = createConfig(docNode,arch,opts,mode) %lib,debug,char,opt,ts)
pc = docNode.createElement('PropertyGroup');
if(mode.debug)
    config = 'Debug';
else
    config = 'Release';
end
pc.setAttribute('Condition',['''$(Configuration)|$(Platform)''==''' config '|' arch '''']);
pc.setAttribute('Label','Configuration');
if(~isempty(opts.configType)),   addElemText(docNode,pc,'ConfigurationType',opts.configType); end
if(mode.debug)
    addElemText(docNode,pc,'UseDebugLibraries','true');
else
    addElemText(docNode,pc,'UseDebugLibraries','false');
end
if(~isempty(opts.charset)),  addElemText(docNode,pc,'CharacterSet',opts.charset); end
if(mode.optimize), addElemText(docNode,pc,'WholeProgramOptimization','true'); end
if(~isempty(opts.toolset)),    addElemText(docNode,pc,'PlatformToolset',opts.toolset); end

%Create Import Section
function pc = createImport(docNode,name,cond,label)
pc = docNode.createElement('Import');
pc.setAttribute('Project',name);
if(~isempty(cond)); pc.setAttribute('Condition',cond); end
if(~isempty(label)); pc.setAttribute('Label',label); end

%Create Property Sheet Section
function pc = createPropSheet(docNode,config,plat)
pc = docNode.createElement('ImportGroup');
pc.setAttribute('Label','PropertySheets');
pc.setAttribute('Condition',['''$(Configuration)|$(Platform)''==''' config '|' plat '''']);
pc.appendChild(createImport(docNode,'$(UserRootDir)\Microsoft.Cpp.$(Platform).user.props','exists(''$(UserRootDir)\Microsoft.Cpp.$(Platform).user.props'')','LocalAppDataPlatform'));

%Write Project Configuration Section
function pc = writeProjConfig(docNode,config,plat)
pc = docNode.createElement('ProjectConfiguration');
pc.setAttribute('Include',[config '|' plat]);
addElemText(docNode,pc,'Configuration',config)
addElemText(docNode,pc,'Platform',plat);

%Write Fortran Configuration Section
function pc = writeIFortConfig(docNode,config,plat,exPP,exInc,ifortver)
pc = docNode.createElement('Configuration');
pc.setAttribute('Name',[config '|' plat]);
if(strcmpi(plat,'x64'))
    pc.setAttribute('OutputDirectory','$(SolutionDir)x64\$(ConfigurationName)');
else
    pc.setAttribute('OutputDirectory','$(SolutionDir)$(ConfigurationName)');
end
pc.setAttribute('ConfigurationType','typeStaticLibrary');
tl = docNode.createElement('Tool');
tl.setAttribute('Name','VFFortranCompilerTool');
tl.setAttribute('SuppressStartupBanner','true');
tl.setAttribute('Preprocess','preprocessYes');
tl.setAttribute('EnableEnhancedInstructionSet','codeArchSSE2');
if(nargin > 3 && ~isempty(exPP))
    tl.setAttribute('PreprocessorDefinitions',concatenatePP(exPP));
end	
if(nargin > 4 && ~isempty(exInc))
    tl.setAttribute('AdditionalIncludeDirectories',concatenatePP(exInc));
end
if(strcmpi(config,'debug'))
    tl.setAttribute('DebugInformationFormat','debugEnabled');
    tl.setAttribute('Optimization','optimizeDisabled');
    tl.setAttribute('RuntimeLibrary','rtMultiThreadedDebugDLL');
    tl.setAttribute('Traceback','true');
    tl.setAttribute('BoundsCheck','true');
    tl.setAttribute('StackFrameCheck','true');
else
    tl.setAttribute('Optimization','optimizeFull');
    tl.setAttribute('RuntimeLibrary','rtMultiThreadedDLL');
end
pc.appendChild(tl);
if(strcmpi(ifortver,'XE15'))
    tl = docNode.createElement('Tool'); %seem to need all these now
    tl.setAttribute('Name','VFLibrarianTool');
    pc.appendChild(tl);
    tl = docNode.createElement('Tool'); 
    tl.setAttribute('Name','VFResourceCompilerTool');
    pc.appendChild(tl);
    tl = docNode.createElement('Tool'); 
    tl.setAttribute('Name','VFMidlTool');
    tl.setAttribute('SuppressStartupBanner','true');
    if(strcmpi(plat,'x64'))
        tl.setAttribute('TargetEnvironment','midlTargetAMD64');
    end    
    pc.appendChild(tl);
    tl = docNode.createElement('Tool'); 
    tl.setAttribute('Name','VFCustomBuildTool');
    pc.appendChild(tl);
    tl = docNode.createElement('Tool'); 
    tl.setAttribute('Name','VFPreLinkEventTool');
    pc.appendChild(tl);
    tl = docNode.createElement('Tool'); 
    tl.setAttribute('Name','VFPreBuildEventTool');
    pc.appendChild(tl);
    tl = docNode.createElement('Tool'); 
    tl.setAttribute('Name','VFPostBuildEventTool');
    pc.appendChild(tl);
else
    if(strcmpi(plat,'x64'))
        tl = docNode.createElement('Tool');
        tl.setAttribute('Name','VFMidlTool');
        tl.setAttribute('SuppressStartupBanner','true');
        tl.setAttribute('TargetEnvironment','midlTargetAMD64');
        pc.appendChild(tl);
    end
end


%Write C++ Debug Detail Section
function pc = writeDebugDetail(docNode,upath,hdr,plat,exPP,openMP,console,linklib,linkpath,deffile,casCpp,excepExtern)
pc = docNode.createElement('ItemDefinitionGroup');
pc.setAttribute('Condition',['''$(Configuration)|$(Platform)''==''Debug|',plat,'''']);
cl = docNode.createElement('ClCompile');
cl.appendChild(docNode.createElement('PrecompiledHeader'));
addElemText(docNode,cl,'WarningLevel','Level3');
addElemText(docNode,cl,'Optimization','Disabled');
if(openMP), addElemText(docNode,cl,'OpenMPSupport','true'); end
if(~isempty(exPP)), exPP = concatenatePP(exPP); end
if(~console), ppType = '_LIB'; else ppType = '_CONSOLE'; end
if(strcmpi(plat,'x64')), pplat = 'WIN32;WIN64'; else pplat = 'WIN32'; end
addElemText(docNode,cl,'PreprocessorDefinitions',[pplat ';_DEBUG;' ppType ';' exPP '%(PreprocessorDefinitions)']);
addElemText(docNode,cl,'AdditionalIncludeDirectories',[includeStr(upath,hdr) '%(AdditionalIncludeDirectories)']);
if(casCpp), addElemText(docNode,cl,'CompileAs','CompileAsCpp'); end
if(excepExtern), addElemText(docNode,cl,'ExceptionHandling','SyncCThrow'); end
pc.appendChild(cl);
lk = docNode.createElement('Link');
if(console)
    addElemText(docNode,lk,'SubSystem','Console');
else
    addElemText(docNode,lk,'SubSystem','Windows');
end
addElemText(docNode,lk,'GenerateDebugInformation','true');
if(~isempty(linklib))
    addElemText(docNode,lk,'AdditionalDependencies',[concatenatePP(linklib) 'kernel32.lib;user32.lib;gdi32.lib;'...
    'winspool.lib;comdlg32.lib;advapi32.lib;shell32.lib;ole32.lib;oleaut32.lib;uuid.lib;odbc32.lib;odbccp32.lib;%(AdditionalDependencies)']);
    if(strcmpi(plat,'x64'))
        libp = concatenatePP(linkpath,'\x64\Debug');
    else
        libp = concatenatePP(linkpath,'\Debug');
    end
    addElemText(docNode,lk,'AdditionalLibraryDirectories',libp);
end
if(~isempty(deffile))
    addElemText(docNode,lk,'ModuleDefinitionFile',deffile);
end
pc.appendChild(lk);

%Write C++ Release Detail Section
function pc = writeReleaseDetail(docNode,upath,hdr,plat,exPP,openMP,console,linklib,linkpath,deffile,casCpp,excepExtern)
pc = docNode.createElement('ItemDefinitionGroup');
pc.setAttribute('Condition',['''$(Configuration)|$(Platform)''==''Release|',plat,'''']);
cl = docNode.createElement('ClCompile');
cl.appendChild(docNode.createElement('PrecompiledHeader'));
addElemText(docNode,cl,'WarningLevel','Level3');
addElemText(docNode,cl,'Optimization','MaxSpeed');
addElemText(docNode,cl,'IntrinsicFunctions','true');
if(~strcmpi(plat,'x64')), addElemText(docNode,cl,'EnableEnhancedInstructionSet','StreamingSIMDExtensions2'); end
addElemText(docNode,cl,'MultiProcessorCompilation','true');
addElemText(docNode,cl,'RuntimeTypeInfo','true');
addElemText(docNode,cl,'FavorSizeOrSpeed','Speed');
addElemText(docNode,cl,'InlineFunctionExpansion','OnlyExplicitInline');
addElemText(docNode,cl,'RuntimeLibrary','MultiThreadedDLL');
if(openMP), addElemText(docNode,cl,'OpenMPSupport','true'); end
if(~isempty(exPP)), exPP = concatenatePP(exPP); end
if(~console), ppType = '_LIB'; else ppType = '_CONSOLE'; end
if(strcmpi(plat,'x64')), pplat = 'WIN32;WIN64'; else pplat = 'WIN32'; end
addElemText(docNode,cl,'PreprocessorDefinitions',[pplat ';NDEBUG;' ppType ';' exPP '%(PreprocessorDefinitions)']);
addElemText(docNode,cl,'AdditionalIncludeDirectories',[includeStr(upath,hdr) '%(AdditionalIncludeDirectories)']);
if(casCpp), addElemText(docNode,cl,'CompileAs','CompileAsCpp'); end
if(excepExtern), addElemText(docNode,cl,'ExceptionHandling','SyncCThrow'); end
pc.appendChild(cl);
lk = docNode.createElement('Link');
if(console)
    addElemText(docNode,lk,'SubSystem','Console');
    addElemText(docNode,lk,'EnableCOMDATFolding','true');
else
    addElemText(docNode,lk,'SubSystem','Windows');
end
addElemText(docNode,lk,'GenerateDebugInformation','false');
addElemText(docNode,lk,'OptimizeReferences','true');
if(~isempty(linklib))
    addElemText(docNode,lk,'AdditionalDependencies',[concatenatePP(linklib) 'kernel32.lib;user32.lib;gdi32.lib;'...
    'winspool.lib;comdlg32.lib;advapi32.lib;shell32.lib;ole32.lib;oleaut32.lib;uuid.lib;odbc32.lib;odbccp32.lib;%(AdditionalDependencies)']);
    if(strcmpi(plat,'x64'))
        libp = concatenatePP(linkpath,'\x64\Release');
    else
        libp = concatenatePP(linkpath,'\Release');
    end
    addElemText(docNode,lk,'AdditionalLibraryDirectories',libp);
end
if(~isempty(deffile))
    addElemText(docNode,lk,'ModuleDefinitionFile',deffile);
end
pc.appendChild(lk);

%Function to get a GUID from .EXE
function str = getProjGUID()
cdir = cd;
cd(GETCD);
dos(['"' GETCD 'genGUID" 1']);
cd(cdir);
fh = fopen([GETCD 'guids.txt']);
str = fgetl(fh);
if(length(str) < 5)
    error('Error reading GUID');
end
fclose(fh);

%Generate an include string from supplied cell array
function str = includeStr(upath,hdr)
str = [];
for i = 1:size(hdr,1) %for each path
    if(iscell(upath))
        p = upath{hdr{i,3}};
    else
        p = upath;
    end
    %Check we are on the same path as our base one
    h = hdr{i,1};
    if(isempty(strfind(h,p))) %empty means not same base path
        attp = getPathSep(h,p);
    else %OK contained with base path
        if(any(strfind(p,':')))
            ind = strfind(p,filesep); %find end file separator
            len2 = ind(end)+1; %drop last folder (we are always back one)
        else
            len2 = 1;
        end    
        attp = hdr{i,1}(len2:end);
    end
    
    str = [str '..' filesep attp ';']; %#ok<AGROW>
end

%Determine how many ../ we need and get path
function str = getPathSep(h,p)
%Run forward until we find path section that doesn't match
for i = 1:min(length(h),length(p))
    if(~strcmp(h(1:i),p(1:i)))
        break;
    end
end
%Go back until last pathsep
for k = i:-1:1
    if(strcmp(h(k),filesep))
        break;
    end
end
i = k+1;
%Now determine how many file seps in front of us
ind = strfind(p(i:end),filesep);
%Determine how many explicit steps back
indb = strfind(p(i:end),['..' filesep]);
str = '';
for j = 1:(length(ind)-length(indb)*2) %check the *2 (or mod of) with libClp
    str = sprintf('%s..%s',str,filesep);
end
%Concatenate with new bit 
str = [str h(i:end)];

%Get cwd of this file
function str = GETCD()
str = which('VS_WriteProj.m');
ind = strfind(str,filesep);
str(ind(end)+1:end) = [];

%Concatentate preprocessors
function str = concatenatePP(pp,extra)
if(isempty(pp))
    str = [];
    return;
end
if(nargin < 2), extra = ''; end
if(iscell(pp))
    str = [];
    for i = 1:length(pp)
        if(isempty(strfind(pp{i},':')))
            str = sprintf('%s%s%s;',str,pp{i},extra);
        else
            str = sprintf('%s%s;',str,pp{i});
        end
    end
elseif(ischar(pp))
    str = pp;
else
    error('Unknown preprocessor format');
end
