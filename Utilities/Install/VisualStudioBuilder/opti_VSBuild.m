function opti_VSBuild(solver,paths,opts)
%Build a Visual Studio Solution/Project for a Selected Solver
%   
%   This function will attempt to generate an OPTI compatible Visual Studio
%   Project of the selected solver, automatically compile it, and then save
%   the generated libraries and header files into the required OPTI folders
%   in order to compile it as a MEX function. The project files will remain 
%   for you to use if you like outside MATLAB.
%
%   opti_VSBuild(solver_name,opts_structure)
%      - solver_name:       CLP, NOMAD, etc
%      - path:              FULL path to the solver source directory
%      - opts_structure:    structure with the following fields:          
%           expath:         FULL paths to extra required source directories (cell array)
%           vsver:          Visual Studio Version (default find the highest version, otherwise 'VS2013', 'VS2010' and 'VS2012' supported)
%           ifortver:       Intel Fortran Version (default find the highest version, otherwise 'XE15', 'XE13SP1', 'XE13', 'XE11SP1' and 'XE11' also supported)
%           ma57            Include MA57 in build {[]} ('HSL' - HSL Fortran source version, 'Matlab' - DLL included with MATLAB) [see opti_MISC_Install]
%           ma27            Include MA27 in build {[]} ('HSL' - HSL Fortran source version) [see opti_MISC_Install]
%           pardiso         Include PARDISO in build {[]} ('Basel' - DLL from Official PARDISO Site, 'MKL' - Version included with Intel MKL)
%           mumps           Include MUMPS in build {false}
%           linloader       Include IPOPT Linear System Loader DLL Code {false}

%Constants (modify to suit your PC)
vs13dir     = 'C:\Program Files (x86)\Microsoft Visual Studio 12.0\Common7\IDE';
vs12dir     = 'C:\Program Files (x86)\Microsoft Visual Studio 11.0\Common7\IDE';
vs10dir     = 'C:\Program Files (x86)\Microsoft Visual Studio 10.0\Common7\IDE';
if15dir     = 'C:\Program Files (x86)\Intel\Composer XE 2015\compiler';
if13SP1dir  = 'C:\Program Files (x86)\Intel\Composer XE 2013 SP1\compiler';
if13dir     = 'C:\Program Files (x86)\Intel\Composer XE 2013\compiler';
if11SP1dir  = 'C:\Program Files (x86)\Intel\Composer XE 2011 SP1\compiler';
if11dir     = 'C:\Program Files (x86)\Intel\Composer XE 2011\compiler';

%Default Options
if(nargin < 3 || isempty(opts))
    opts.vsver = [];
    opts.ifortver = [];
    opts.pardiso = [];
    opts.ma57 = [];  
    opts.ma27 = [];
    opts.linloader = false;
end

%Special paths
if(nargin < 2)
    switch(lower(solver))
        case 'mwma57'
            paths = {[cd '\Solvers\Source\opti']};
        otherwise
            error('You must supply the FULL path to the solver as the second argument!');
    end
end

%Extract options
if(~iscell(paths)), paths = {paths}; end
if(isfield(opts,'expaths') && ~isempty(opts.expaths))
    if(~iscell(paths))
        paths = [paths {opts.expaths}];
    else
        paths = [paths opts.expaths];
    end   
end
if(~isfield(opts,'pardiso')), PARDISO = []; else PARDISO = opts.pardiso; end
if(~isfield(opts,'ma57')), MA57 = []; else MA57 = opts.ma57; end
if(~isfield(opts,'ma27')), MA27 = []; else MA27 = opts.ma27; end
if(~isfield(opts,'mumps') || isempty(opts.mumps)), MUMPS = []; else MUMPS = opts.mumps; end
if(~isfield(opts,'linloader') || isempty(opts.linloader)), LINLOADER = false; else LINLOADER = opts.linloader; end

%Visual Studio Setup
if(isfield(opts,'vsver') && ~isempty(opts.vsver))
    vsver = opts.vsver;
else %find suitable version
    if(exist([vs13dir '\devenv.exe'],'file'))
        vsver = 'VS2013';
    elseif(exist([vs12dir '\devenv.exe'],'file'))
        vsver = 'VS2012';
    elseif(exist([vs10dir '\devenv.exe'],'file'))
        vsver = 'VS2010';
    else
        error('Cannot find a suitable version of Visual Studio!');
    end
end
%Intel Fortran Setup
if(isfield(opts,'ifortver') && ~isempty(opts.ifortver))
    ifortver = opts.ifortver;
else %find suitable version
    if(exist(if15dir,'dir'))
        ifortver = 'XE15';
    elseif(exist(if13dir,'dir') || exist(if13SP1dir,'dir'))
        ifortver = 'XE13';
    elseif(exist(if11dir,'dir') || exist(if11SP1dir,'dir'))
        ifortver = 'XE11';
    else
        ifortver = [];
    end
end
if(~isempty(ifortver)) %ensure compatible with visual studio builder
    switch(lower(ifortver))
        case {'11.2','xe15'}
            ifortver = 'XE15';
        case {'11.1','11.0','11','xe13','xe13sp1'}
            ifortver = 'XE13';
        case {'10.3sp1','10.3 sp1','xe11','xe11sp1'}
            ifortver = 'XE11';
    end
end
%Current working directory
cdir = cd;

%Find Visual Studio devenv.exe 
switch(lower(vsver))
    case {'vs2013','2013'}
        if(~exist([vs13dir '\devenv.exe'],'file'))
            error('Could not find Visual Studio 2013 Install - Looked in ''%s''',vs13dir);
        end
        vsdir = vs13dir;        
    case {'vs2012','2012'}
        if(~exist([vs12dir '\devenv.exe'],'file'))
            error('Could not find Visual Studio 2012 Install - Looked in ''%s''',vs12dir);
        end
        vsdir = vs12dir;        
    case {'vs2010','2010'}
        if(~exist([vs10dir '\devenv.exe'],'file'))
            error('Could not find Visual Studio 2010 Install - Looked in ''%s''',vs10dir);
        end
        vsdir = vs10dir;        
    otherwise
        error('Unknown Visual Studio Version - Use VS2010, VS2012 or VS2013');
end
%Misc Defines
bthdr = [cdir '/Solvers/Source/Include/BuildTools']; %Build Tools Headers for Coin Solvers

%Check at least first path can be found
if(~exist(paths{1},'dir'))
    error('Cannot find directory ''%s''!',paths{1});
end

switch(lower(solver))
    case 'clp'
        clppath = paths{1};
        if(length(paths) >= 2)
            glpkpath = paths{2};
        else
            glpkpath = [];
        end
        %Copy Header Files
        if(~exist([cdir '/Solvers/Source/Include/Clp'],'dir'))
            mkdir([cdir '/Solvers/Source/Include/Clp']);
            mkdir([cdir '/Solvers/Source/Include/Coin']);
            mkdir([cdir '/Solvers/Source/Include/Osi']);
        end
        %Clp Headers
        [~,hdrs] = VS_BuildFileList([clppath '/src']);
        copyHeaders(hdrs,[cdir '/Solvers/Source/Include/Clp/']);
        %Coin Headers
        [~,hdrs] = VS_BuildFileList([clppath '/../CoinUtils/src']);
        copyHeaders(hdrs,[cdir '/Solvers/Source/Include/Coin/']);
        %OSI Headers
        [~,hdrs] = VS_BuildFileList([clppath '/../Osi/src/Osi']);
        copyHeaders(hdrs,[cdir '/Solvers/Source/Include/Osi/']);
        n = 1;
        %CLP
        sdir = [clppath '\src'];
        hdrs = {[clppath '\..\CoinUtils\src'], [clppath '\..\Osi\src'], bthdr}; 
        name = 'libclp';
        opts = [];
        opts.exPP = {'CLP_BUILD','_CRT_SECURE_NO_WARNINGS'};
        opts.exclude = {'MyEventHandler.cpp','MyMessageHandler.cpp','unitTest.cpp','CbcOrClpParam.cpp',...
                        'ClpCholeskyMumps.cpp','ClpCholeskyUfl.cpp','ClpCholeskyWssmp.cpp','ClpCholeskyWssmpKKT.cpp','ClpMain.cpp'};
        opts.exFilter = {'Abc*','CoinAbc*'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %CLP with Aboca (requires Intel C++ Compiler)
        sdir = [clppath '\src'];
        hdrs = {[clppath '\..\CoinUtils\src'], [clppath '\..\Osi\src'], bthdr}; 
        name = 'libclpabc';
        opts = [];
        opts.exPP = {'CLP_BUILD','_CRT_SECURE_NO_WARNINGS','CLP_HAS_ABC=4','__BYTE_ORDER=__LITTLE_ENDIAN','INTEL_COMPILER'};
        opts.exclude = {'MyEventHandler.cpp','MyMessageHandler.cpp','unitTest.cpp','CbcOrClpParam.cpp',...
                        'ClpCholeskyMumps.cpp','ClpCholeskyUfl.cpp','ClpCholeskyWssmp.cpp','ClpCholeskyWssmpKKT.cpp','ClpMain.cpp'};
        opts.toolset = 'IntelC++';
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        % CoinUtils
        sdir = [clppath '\..\CoinUtils\src'];
        hdrs = {bthdr};
        name = 'libcoinutils';
        opts = [];
        opts.exPP = {'COINUTILS_BUILD','_CRT_SECURE_NO_WARNINGS'};          
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        % CoinUtils with GLPK (for GMPL CoinR reading)
        if(~isempty(glpkpath))
            sdir = [clppath '\..\CoinUtils\src'];
            hdrs = {[glpkpath '\src'], bthdr};
            name = 'libcoinutilsgmpl';
            opts = [];
            opts.exPP = {'COIN_HAS_GLPK','COINUTILS_BUILD','_CRT_SECURE_NO_WARNINGS'};          
            VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        end
        % Osi
        sdir = [clppath '\..\Osi\src'];
        hdrs = {[clppath '\..\CoinUtils\src'] bthdr};
        name = 'libosi';
        opts = [];
        opts.exPP = {'OSI_BUILD','_CRT_SECURE_NO_WARNINGS'};   
        opts.exFolder = {'OsiCpx','OsiGlpk','OsiGrb','OsiMsk','OsiSpx','OsiXpr','OsiCommonTest'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        % Clp Solver (optional)
        sdir = [clppath '\src'];
        hdrs = {[clppath '\..\CoinUtils\src'], [clppath '\..\Osi\src'], bthdr};
        name = 'Clp';
        opts = [];
        opts.exPP = {'CLP_BUILD','_CRT_SECURE_NO_WARNINGS'};
        opts.include = {'ClpMain.cpp','CbcOrClpParam.cpp','unitTest.cpp','MyEventHandler.cpp','MyMessageHandler.cpp'};
        opts.console = true;
        opts.linklib = {'libclp.lib','libcoinutils.lib','libosi.lib'};
        opts.linkpath = {'..\libclp'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'libcoinutils','libosi','libclp'};  %,'libclpabc'
        comps = {'vc','vc','vc'};
        if(~isempty(glpkpath))
            projs = [projs 'libcoinutilsgmpl'];
            comps = [comps 'vc'];
        end
        
    case 'cbc'
        cbcpath = paths{1};
        %Copy Header Files
        if(~exist([cdir '/Solvers/Source/Include/Cbc'],'dir'))
            mkdir([cdir '/Solvers/Source/Include/Cbc']);
            mkdir([cdir '/Solvers/Source/Include/Cgl']);
        end
        %Cbc Headers
        [~,hdrs] = VS_BuildFileList([cbcpath '/src']);
        copyHeaders(hdrs,[cdir '/Solvers/Source/Include/Cbc/']);
        %Cgl Headers
        [~,hdrs] = VS_BuildFileList([cbcpath '/../Cgl/src']);
        copyHeaders(hdrs,[cdir '/Solvers/Source/Include/Cgl/']);
        n = 1;
        % CBC
        sdir = [cbcpath '\src'];
        hdrs = {[cbcpath '\..\CoinUtils\src'], [cbcpath '\..\Cgl\src'], [cbcpath '\..\Clp\src'], [cbcpath '\..\Osi\src'], bthdr};
        name = 'libcbc';
        opts = [];
        opts.exPP = {'CBC_BUILD','USE_CBCCONFIG','COIN_NO_TEST_DUPLICATE','_CRT_SECURE_NO_WARNINGS'};
        opts.exclude = {'unitTest.cpp','CoinSolve.cpp','CbcGeneric.cpp','CbcBranchBase.cpp'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        % CGL
        sdir = [cbcpath '\..\Cgl\src'];
        hdrs = {[cbcpath '\..\CoinUtils\src'], [cbcpath '\..\Clp\src\'], [cbcpath '\..\Osi\src'], bthdr};
        name = 'libcgl';
        opts = [];
        opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};          
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        % OPTIONAL PARTS
        %CLP
        sdir = [cbcpath '\..\Clp\src'];
        hdrs = {[cbcpath '\..\CoinUtils\src'], [cbcpath '\..\Osi\src'], bthdr};
        name = 'libclp';
        opts = [];
        opts.exPP = {'CLP_BUILD','_CRT_SECURE_NO_WARNINGS'};
        opts.exclude = {'MyEventHandler.cpp','MyMessageHandler.cpp','unitTest.cpp','CbcOrClpParam.cpp',...
                        'ClpCholeskyMumps.cpp','ClpCholeskyUfl.cpp','ClpCholeskyWssmp.cpp','ClpCholeskyWssmpKKT.cpp','ClpMain.cpp'};
        opts.exFilter = {'Abc*','CoinAbc*'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %CLP with Aboca (requires Intel C++ Compiler)
        sdir = [cbcpath '\src'];
        hdrs = {[cbcpath '\..\CoinUtils\src'], [cbcpath '\..\Osi\src'], bthdr};
        name = 'libclpabc';
        opts = [];
        opts.exPP = {'CLP_BUILD','_CRT_SECURE_NO_WARNINGS','CLP_HAS_ABC=4','__BYTE_ORDER=__LITTLE_ENDIAN','INTEL_COMPILER'};
        opts.exclude = {'MyEventHandler.cpp','MyMessageHandler.cpp','unitTest.cpp','CbcOrClpParam.cpp',...
                        'ClpCholeskyMumps.cpp','ClpCholeskyUfl.cpp','ClpCholeskyWssmp.cpp','ClpCholeskyWssmpKKT.cpp','ClpMain.cpp'};
        opts.toolset = intc;
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        % CoinUtils
        sdir = [cbcpath '\..\CoinUtils\src'];
        hdrs = {bthdr};
        name = 'libcoinutils';
        opts = [];
        opts.exPP = {'COINUTILS_BUILD','_CRT_SECURE_NO_WARNINGS'};          
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        % Osi
        sdir = [cbcpath '\..\Osi\src'];
        hdrs = {[cbcpath '\..\CoinUtils\src'] bthdr};
        name = 'libosi';
        opts = [];
        opts.exPP = {'OSI_BUILD','_CRT_SECURE_NO_WARNINGS'};   
        opts.exFolder = {'OsiCpx','OsiGlpk','OsiGrb','OsiMsk','OsiSpx','OsiXpr','OsiCommonTest'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        % CBC Solver (optional)
        sdir = [cbcpath '\src'];
        hdrs = {[cbcpath '\..\CoinUtils\src'], [cbcpath '\..\Cgl\src'], [cbcpath '\..\Clp\src'], [cbcpath '\..\Osi\src'], bthdr};
        name = 'Cbc';
        opts = [];
        opts.exPP = {'CBC_BUILD','USE_CBCCONFIG','COIN_NO_TEST_DUPLICATE','_CRT_SECURE_NO_WARNINGS'};
        opts.include = {'CoinSolve.cpp'};
        opts.console = true;
        opts.linklib = {'libcbc.lib','libclp.lib','libcgl.lib','libcoinutils.lib','libosi.lib'};
        opts.linkpath = {'..\libcbc'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'libcbc','libcgl'};  
        comps = {'vc','vc'};
        
    case 'metis'
        metispath = paths{1};
        %Comment strings.h
        str = fileread([metispath '/Lib/metis.h']);
        if(isempty(strfind(str,'//#include <strings.h>')))
            str = regexprep(str,'#include <strings.h>','//#include <strings.h>');
            fp = fopen([metispath '/Lib/metis.h'],'w');
            fprintf(fp,'%s',str); fclose(fp);
        end
        %   	- Open proto.h and uncomment line 435 (void GKfree...)
        %       - (METIS >= v5 only) in metis.h define REALTYPEWIDTH 64
        %       - (METIS >= v5 only) in metis.h define USE_GKREGEX
        %Copy Header Files
        if(~exist([cdir '/Solvers/Source/Include/metis'],'dir'))            
            mkdir([cdir '/Solvers/Source/Include/metis']);
        end
        %Metis Headers
        [~,hdrs] = VS_BuildFileList([metispath '/Lib']);
        copyHeaders(hdrs,[cdir '/Solvers/Source/Include/metis/']);
        % METIS
        n = 1;
        sdir = [metispath '\Lib'];
        name = 'libmetis';
        opts = [];
        opts.exPP = {'_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE','__STDC__','__VC__'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'libmetis'};  
        comps = {'vc'};
        
    case 'mumps'
        mumpspath = paths{1};
        if(length(paths) > 1) %assume metis path 2nd
            opti_VSBuild('metis',paths{2}); %build METIS first
        end
        %Copy Header Files
        if(~exist([cdir '/Solvers/Source/Include/Mumps'],'dir'))
            mkdir([cdir '/Solvers/Source/Include/Mumps']);
        end
        %Mumps Headers
        [~,hdrs] = VS_BuildFileList([mumpspath '/include']);
        copyHeaders(hdrs,[cdir '/Solvers/Source/Include/Mumps/']);       
        n = 1;
        % DMUMPS_C
        sdir = [mumpspath '\src'];
        hdrs = {[mumpspath '\include'],[mumpspath '\PORD\include']};
        name = 'libdmumps_c';
        opts = [];
        opts.exPP = {'_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE','MUMPS_ARITH=MUMPS_ARITH_d','pord','metis'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        % DMUMPS_F
        sdir = [mumpspath '\src'];
        hdrs = {[mumpspath '\include'],[mumpspath '\libseq']};
        name = 'libdmumps_f';
        opts = [];
        opts.cpp = false; opts.ifortver = ifortver;
        opts.exPP = {'pord','metis'};
        opts.exFilter = {'cmumps*','smumps*','zmumps*'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        % ZMUMPS_C 
        sdir = [mumpspath '\src'];
        hdrs = {[mumpspath '\include'],[mumpspath '\PORD\include']};
        name = 'libzmumps_c';
        opts = [];
        opts.exPP = {'_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE','MUMPS_ARITH=MUMPS_ARITH_z','pord','metis'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        % ZMUMPS_F 
        sdir = [mumpspath '\src'];
        hdrs = {[mumpspath '\include'],[mumpspath '\libseq']};
        name = 'libzmumps_f';
        opts = [];
        opts.cpp = false; opts.ifortver = ifortver;
        opts.exPP = {'pord','metis'};
        opts.exFilter = {'cmumps*','smumps*','dmumps*'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        % LIBSEQ_C
        sdir = [mumpspath '\libseq'];
        name = 'libseq_c';
        opts = [];
        opts.exPP = {'_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        % LIBSEQ_F
        sdir = [mumpspath '\libseq'];
        name = 'libseq_f';
        opts = [];       
        opts.cpp = false; opts.ifortver = ifortver;
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        % PORD
        sdir = [mumpspath '\PORD\lib'];
        hdrs = [mumpspath '\PORD\include'];
        name = 'libpord';
        opts = [];
        opts.exPP = {'_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'libdmumps_c','libzmumps_c','libseq_c','libpord','libdmumps_f','libzmumps_f','libseq_f'};  
        comps = {'vc','vc','vc','if','if','if','if'};
        
    case 'ipopt'
        ipoptpath = paths{1};
        if(length(paths) > 1)
            mumpspath = paths{2};
        else
            mumpspath = [];
        end
        %Comment hsl & asl defines
        str = fileread([ipoptpath '/src/Common/config_default.h']);
        if(isempty(strfind(str,'//#define COIN_HAS_ASL 1')))
            str = regexprep(str,'#define COIN_HAS_ASL 1','//#define COIN_HAS_ASL 1');
            str = regexprep(str,'#define COIN_HAS_HSL 1','//#define COIN_HAS_HSL 1');
            fp = fopen([ipoptpath '/src/Common/config_default.h'],'w');
            fprintf(fp,'%s',str); fclose(fp);
        end
        %Copy Header Files
        if(~exist([cdir '/Solvers/Source/Include/Ipopt'],'dir'))
            mkdir([cdir '/Solvers/Source/Include/Ipopt']);
        end
        %Ipopt Headers
        [~,hdrs] = VS_BuildFileList([ipoptpath '/src']);
        copyHeaders(hdrs,[cdir '/Solvers/Source/Include/Ipopt/']);  
        n=1;
        %Build IPOPT
        sdir = [ipoptpath '\src'];
        hdrs = {bthdr};
        name = 'libipopt';
        opts = [];
        opts.exFolder = {}; opts.exclude = [];
        opts.exPP = {'IPOPT_BUILD','_CRT_SECURE_NO_WARNINGS'};
        if(~LINLOADER) %can't have both dynamic and static linking (for now)
            if(~isempty(PARDISO))
                switch(lower(PARDISO))
                    case 'basel'
                        opts.exPP = [opts.exPP,'HAVE_PARDISO','HAVE_PARDISO_PARALLEL'];
                    case 'mkl'
                        opts.exPP = [opts.exPP,'HAVE_PARDISO','HAVE_PARDISO_OLDINTERFACE','HAVE_PARDISO_MKL'];
                end
            end
            if(~isempty(MA57))
                switch(lower(MA57))
                    case 'hsl'
                        opts.exPP = [opts.exPP,'COINHSL_HAS_MA57'];
                    case 'matlab'
                        checkLibMWMA57;
                        opts.exPP = [opts.exPP,'FUNNY_MA57_FINT','COINHSL_HAS_MA57'];
                end
            end
            if(~isempty(MA27))
                switch(lower(MA27))
                    case 'hsl'
                        opts.exPP = [opts.exPP,'COINHSL_HAS_MA27'];
                end
            end
            if(~isempty(mumpspath) && MUMPS)
                checkLibMUMPS(paths(2:end));
                hdrs = [hdrs [mumpspath '\include'], [mumpspath '\libseq']];
                opts.exPP = [opts.exPP,'COIN_HAS_METIS','COIN_HAS_MUMPS'];
            else
                opts.exclude = {'IpMumpsSolverInterface.cpp'};
            end
            opts.exFolder = {'contrib\LinearSolverLoader'};
        else
            hdrs = [hdrs [ipoptpath '\..\ThirdParty\HSL']];
            opts.exPP = [opts.exPP,'HAVE_LINEARSOLVERLOADER','HAVE_WINDOWS_H','SHAREDLIBEXT="dll"','COINHSL_HSL2013'];
            %post = [post ' -L' libdir ' -llibipoptdyn -llibma28part -DHAVE_LINEARSOLVERLOADER'];
            opts.charset = 'multibyte';    
            opts.outname = 'libipoptdyn';
        end
        opts.exFolder = [opts.exFolder,'Apps','Algorithm\Inexact']; 
        if(~isempty(opts.exclude)) %assume cell
            opts.exclude = [opts.exclude 'AmplTNLP.cpp'];
        else
            opts.exclude = 'AmplTNLP.cpp';
        end
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'libipopt'};
        comps = {'vc'};
        
    case 'bonmin'
        bminpath = paths{1};
        ipoptpath = paths{2};
        mumpspath = paths{3};
        cbcpath = paths{4};
        clppath = paths{5};
        %Comment asl define
        str = fileread([bminpath '/src/Interfaces/config_default.h']);
        if(isempty(strfind(str,'//#define COIN_HAS_ASL 1')))
            str = regexprep(str,'#define COIN_HAS_ASL 1','//#define COIN_HAS_ASL 1');
            fp = fopen([bminpath '/src/Interfaces/config_default.h'],'w');
            fprintf(fp,'%s',str); fclose(fp);
        end
        %Copy Header Files
        if(~exist([cdir '/Solvers/Source/Include/Bonmin'],'dir'))
            mkdir([cdir '/Solvers/Source/Include/Bonmin']);
        end
        %Bonmin Headers
        [~,hdrs] = VS_BuildFileList([bminpath '/src']);
        copyHeaders(hdrs,[cdir '/Solvers/Source/Include/Bonmin/']);  
        n = 1;
        %BONMIN
        sdir = [bminpath '\src'];
        hdrs = {[cbcpath '\..\Cgl\src'],[cbcpath '\src'],[clppath '\src'],[clppath '\..\CoinUtils\src'],[ipoptpath '\src'],[clppath '\..\Osi\src'], bthdr};
        name = 'libbonmin';
        opts = [];
        opts.exFolder = {};
        opts.exPP = {'BONMIN_BUILD','IPOPT_BUILD','_CRT_SECURE_NO_WARNINGS'};
        if(~haveLinearSolverLoader) %can't have both dynamic and static linking (for now)
            switch(havePardiso)
                case 'BASEL'
                    opts.exPP = [opts.exPP,'HAVE_PARDISO','HAVE_PARDISO_PARALLEL'];
                case 'MKL'
                    opts.exPP = [opts.exPP,'HAVE_PARDISO','HAVE_PARDISO_OLDINTERFACE','HAVE_PARDISO_MKL'];
            end
            if(haveMA57)
                opts.exPP = [opts.exPP,'FUNNY_MA57_FINT','COINHSL_HAS_MA57'];
            end
            if(~isempty(mumpspath))
                hdrs = [hdrs [mumpspath '\include'], [mumpspath '\libseq']];
                opts.exPP = [opts.exPP,'COIN_HAS_METIS','COIN_HAS_MUMPS'];
            end
            opts.exFolder = [opts.exFolder,'contrib\LinearSolverLoader'];
        else
            hdrs = [hdrs [ipoptpath '\..\ThirdParty\HSL']];
            opts.exPP = [opts.exPP,'HAVE_LINEARSOLVERLOADER','HAVE_WINDOWS_H','SHAREDLIBEXT="dll"','COINHSL_HSL2013'];
            opts.charset = 'multibyte';    
            opts.outname = 'libipoptdyn';
        end
        opts.exclude = {'BonCurvatureEstimator.cpp','BonCurvBranchingSolver.cpp'};
        opts.exFolder = {'Apps','Interfaces\Ampl','Interfaces\Filter','Algorithms\Ampl'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;        
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'libbonmin'};  
        comps = {'vc'};
        
    case 'dsdp'
        dsdppath = paths{1};
        %Copy Header Files
        if(~exist([cdir '/Solvers/Source/Include/Dsdp'],'dir'))
            mkdir([cdir '/Solvers/Source/Include/Dsdp']);
        end
        %DSDP Headers
        [~,hdrs] = VS_BuildFileList([dsdppath '/include']);
        copyHeaders(hdrs,[cdir '/Solvers/Source/Include/Dsdp/']);  
        n = 1;
        sdir = [dsdppath '\src'];
        hdrs = [dsdppath '\include'];
        name = 'libdsdp';
        opts = [];
        opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;        
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'libdsdp'};  
        comps = {'vc'};

    case 'filtersd'
        fsdpath = paths{1};
        %Copy over OPTI modifications
        unzip([cd '/Solvers/Source/filterSD/filterSD_JCEdit.zip'],[cd '/Solvers/Source/filterSD/jc']);
        pause(0.1); rehash;
        copyfile([cd '/Solvers/Source/filterSD/jc/filterSD.f'],[fsdpath '\filterSD.f'],'f');
        copyfile([cd '/Solvers/Source/filterSD/jc/glcpd.f'],[fsdpath '\glcpd.f'],'f');
        rmdir([cd '/Solvers/Source/filterSD/jc'],'s');
        n = 1;
        %FilterSD (dense)
        sdir = fsdpath;
        name = 'libfilterSD';
        opts = [];       
        opts.cpp = false; opts.ifortver = ifortver;
        opts.include = {'checkd.f','filterSD.f','glcpd.f','l1sold.f','denseL.f','denseA.f','util.f'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        % FILTERSD (sparse)
        sdir = fsdpath;
        name = 'libfilterSDsp';
        opts = [];       
        opts.cpp = false; opts.ifortver = ifortver;
        opts.include = {'checkd.f','filterSD.f','glcpd.f','l1sold.f','schurQR.f','sparseA.f','util.f'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'libfilterSD','libfilterSDsp'};  
        comps = {'if','if'};
        
    case 'glpk'
        glpkpath = paths{1};
        %Copy Header Files
        if(~exist([cdir '/Solvers/Source/Include/Glpk'],'dir'))
            mkdir([cdir '/Solvers/Source/Include/Glpk']);
        end
        %GLPK Headers
        [~,hdrs] = VS_BuildFileList([glpkpath '/src']);
        copyHeaders(hdrs,[cdir '/Solvers/Source/Include/Glpk/']); 
        n=1;
        sdir = [glpkpath '\src'];
        hdrs = [];
        name = 'libglpk';
        opts = [];
        opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'libglpk'};  
        comps = {'vc'};
        
    case 'minpack'
        n = 1;
        sdir = paths{1};
        name = 'libminpack';
        opts = [];       
        opts.cpp = false; opts.ifortver = ifortver;
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'libminpack'};  
        comps = {'if'};
        
    case 'lbfgsb'
        n = 1;
        sdir = paths{1};
        name = 'liblbfgsb';
        opts = [];       
        opts.cpp = false; opts.ifortver = ifortver;
        opts.include = {'lbfgsb.f','timer.f','linpack.f'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'liblbfgsb'};  
        comps = {'if'};
        
    case 'levmar'
        levpath = paths{1};
        %Copy Header Files
        if(~exist([cdir '/Solvers/Source/Include/Levmar'],'dir'))
            mkdir([cdir '/Solvers/Source/Include/Levmar']);
        end
        %GLPK Headers
        [~,hdrs] = VS_BuildFileList(levpath);
        copyHeaders(hdrs,[cdir '/Solvers/Source/Include/Levmar/']); 
        n = 1;
        sdir = levpath;
        name = 'liblevmar';
        opts = [];
        opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};
        opts.exclude = {'expfit.c','Axb_core.c','lm_core.c','lmbc_core.c','lmblec_core.c','lmbleic_core.c','lmlec_core.c','misc_core.c'};
        opts.exFolder = {'matlab'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'liblevmar'};  
        comps = {'vc'};
        
    case 'lpsolve'
        lpspath = paths{1};
        %Copy Header Files
        if(~exist([cdir '/Solvers/Source/Include/Lpsolve'],'dir'))
            mkdir([cdir '/Solvers/Source/Include/Lpsolve']);
        end
        %LPSOLVE Headers
        [~,hdrs] = VS_BuildFileList(lpspath,{'bfp','colamd','demo','shared','lp_solve','lpsolve55'});
        copyHeaders(hdrs,[cdir '/Solvers/Source/Include/Lpsolve/']);         
        %Copy lusol.c + .h
        copyfile([lpspath '\bfp\bfp_LUSOL\LUSOL\lusol.c'],[lpspath '\lusol.c'],'f');
        copyfile([lpspath '\bfp\bfp_LUSOL\LUSOL\lusol.h'],[lpspath '\lusol.h'],'f');
        n = 1;
        sdir = lpspath;
        hdrs = [lpspath '\bfp\bfp_LUSOL\LUSOL'];
        name = 'liblpsolve';
        opts = [];
        opts.exPP = {'PARSER_LP','INVERSE_ACTIVE=INVERSE_LUSOL','RoleIsExternalInvEngine',...
                     '_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE','LoadableBlasLib=0'};
        opts.exclude = {'lp_solveDLL.c'};
        opts.exFilter = {'lp_BFP*'};
        opts.exFolder = {'lpsolve55','lp_solve','demo','bfp\bfp_LUSOL\LUSOL'}; 
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'liblpsolve'};  
        comps = {'vc'};
        
    case 'm1qn3'
        m1qpath = paths{1};
        n = 1;
        sdir = [m1qpath '\src'];
        name = 'libm1qn3';
        opts = [];       
        opts.cpp = false; opts.ifortver = ifortver;        
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'libm1qn3'};  
        comps = {'if'};
        
    case 'nl2sol'
        portpath = paths{1};
        %Fix machine constants
        str = fileread([portpath '/Mach/d1mach.f']);
        str = strrep(str,sprintf('\n      DATA SMALL(1),SMALL(2) /    1048576,          0 /'),sprintf('\nC      DATA SMALL(1),SMALL(2) /    1048576,          0 /'));
        str = strrep(str,sprintf('\n      DATA LARGE(1),LARGE(2) / 2146435071,         -1 /'),sprintf('\nC      DATA LARGE(1),LARGE(2) / 2146435071,         -1 /'));
        str = strrep(str,sprintf('\n      DATA RIGHT(1),RIGHT(2) / 1017118720,          0 /'),sprintf('\nC      DATA RIGHT(1),RIGHT(2) / 1017118720,          0 /'));
        str = strrep(str,sprintf('\n      DATA DIVER(1),DIVER(2) / 1018167296,          0 /'),sprintf('\nC      DATA DIVER(1),DIVER(2) / 1018167296,          0 /'));
        str = strrep(str,sprintf('\n      DATA LOG10(1),LOG10(2) / 1070810131, 1352628735 /, SC/987/'),sprintf('\nC      DATA LOG10(1),LOG10(2) / 1070810131, 1352628735 /, SC/987/'));
        str = strrep(str,'C      DATA SMALL(1),SMALL(2) /          0,    1048576 /','      DATA SMALL(1),SMALL(2) /          0,    1048576 /');
        str = strrep(str,'C      DATA LARGE(1),LARGE(2) /         -1, 2146435071 /','      DATA LARGE(1),LARGE(2) /         -1, 2146435071 /');
        str = strrep(str,'C      DATA RIGHT(1),RIGHT(2) /          0, 1017118720 /','      DATA RIGHT(1),RIGHT(2) /          0, 1017118720 /');
        str = strrep(str,'C      DATA DIVER(1),DIVER(2) /          0, 1018167296 /','      DATA DIVER(1),DIVER(2) /          0, 1018167296 /');
        str = strrep(str,'C      DATA LOG10(1),LOG10(2) / 1352628735, 1070810131 /, SC/987/','      DATA LOG10(1),LOG10(2) / 1352628735, 1070810131 /, SC/987/');
        fp = fopen([portpath '/Mach/d1mach.f'],'w');
        fprintf(fp,'%s',str); fclose(fp);
        n=1;
        sdir = portpath;
        name = 'libnl2sol';
        opts = [];       
        opts.cpp = false; opts.ifortver = ifortver;
        opts.include = {'dn2f.f','dn2g.f','dn2fb.f','dn2gb.f','dn2rdp.f','dv7scp.f','divset.f','Mach\d1mach.f','Mach\r1mach.f',...
                        'Mach\i1mach.f','drn2g.f','drn2gb.f','dd7tpr.f','dd7upd.f','dg7itb.f','ditsum.f','dl7vml.f',...
                        'dq7apl.f','dq7rad.f','dr7tvm.f','dv7cpy.f','dv2nrm.f','dc7vfn.f','dg7lit.f','dn2cvp.f',...
                        'dn2lrd.f','da7sst.f','df7dhb.f','df7hes.f','dg7qsb.f','dg7qts.f','dl7itv.f','dl7ivm.f',...
                        'dl7msb.f','dl7mst.f','dl7nvr.f','dl7sqr.f','dl7srt.f','dl7svn.f','dl7svx.f','dl7tsq.f',...
                        'dl7tvm.f','do7prd.f','dparck.f','dq7rsh.f','dr7mdc.f','drldst.f','ds7dmp.f','ds7ipr.f',...
                        'ds7lup.f','ds7lvm.f','dv2axy.f','dv7dfl.f','dv7ipr.f','dv7scl.f','dv7vmp.f','i7copy.f',...
                        'i7mdcn.f','i7pnvr.f','i7shft.f','stopx.f','dh2rfa.f','dh2rfg.f','ds7bqn.f','fdump.f',...
                        'seterr.f','dd7mlp.f','dv7shf.f','e9rint.f','eprint.f','i8save.f','s88fmt.f','sdump.f',...
                        'stkdmp.f','a9rntc.f','a9rntd.f','a9rnti.f','a9rntl.f','a9rntr.f','frmatd.f','frmati.f',...
                        'frmatr.f','i0tk00.f','u9dmp.f','i10wid.f','iceil.f','iflr.f'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'libnl2sol'};  
        comps = {'if'};
        
    case 'nlopt'
        nloptpath = paths{1};
        %Copy Header Files
        if(~exist([cdir '/Solvers/Source/Include/Nlopt'],'dir'))
            mkdir([cdir '/Solvers/Source/Include/Nlopt']);
        end
        %Nlopt Headers
        [~,hdrs] = VS_BuildFileList([nloptpath '/api']);
        copyHeaders(hdrs,[cdir '/Solvers/Source/Include/Nlopt/']);  
        %Copy config.h
        if(strcmpi(vsver,'VS2013'))
            copyfile([cd '/Solvers/Source/nlopt/config.h'],[nloptpath '\api\config.h'],'f');
        else
            copyfile([cd '/Solvers/Source/nlopt/config12.h'],[nloptpath '\api\config.h'],'f');
        end
        %VS2013 no-vector fix for cobyla.c (see https://connect.microsoft.com/VisualStudio/feedback/details/1028781/c1001-on-release-build)
        if(strcmpi(vsver,'VS2013'))
            str = fileread([nloptpath '/cobyla/cobyla.c']);
            str = regexprep(str,'i__1 = nact;\n  for (k = 1;','i__1 = nact;\n  #pragma loop(no_vector)\n  for (k = 1;');
            fp = fopen([nloptpath '/cobyla/cobyla.c'],'w');
            fprintf(fp,'%s',str); fclose(fp);
        end
        n=1;
        sdir = nloptpath;
        name = 'libnlopt';
        opts = [];
        opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};
        opts.exclude = {'testfuncs.c','tst.cc','tstc.c','testros.cc','prog.cc','redblack_test.c','DIRparallel.c'};
        opts.exFolder = {'octave','swig','test'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'libnlopt'};  
        comps = {'vc'};
        
    case 'nomad'
        nmdpath = paths{1};
        %Copy Header Files
        if(~exist([cdir '/Solvers/Source/Include/Nomad'],'dir'))
            mkdir([cdir '/Solvers/Source/Include/Nomad']);
        end
        %Nlopt Headers
        [~,hdrs] = VS_BuildFileList(nmdpath);
        copyHeaders(hdrs,[cdir '/Solvers/Source/Include/Nomad/']);
        %Increase max dim
        str = fileread([nmdpath '/defines.hpp']);
        if(isempty(strfind(str,'const int MAX_DIMENSION = 5000;')))
            str = regexprep(str,'const int MAX_DIMENSION = 1000;','const int MAX_DIMENSION = 5000;');
            fp = fopen([nmdpath '/defines.hpp'],'w');
            fprintf(fp,'%s',str); fclose(fp);
        end
        n=1;
        sdir = nmdpath;
        name = 'libnomad';
        opts = [];
        opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};
        opts.exclude = 'nomad.cpp';
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'libnomad'};  
        comps = {'vc'};
        
    case 'ooqp'       
        ooqppath = paths{1};
        %Copy Header Files
        if(~exist([cdir '/Solvers/Source/Include/Ooqp'],'dir'))
            mkdir([cdir '/Solvers/Source/Include/Ooqp']);
        end
        %OOQP Headers
        [~,hdrs] = VS_BuildFileList(ooqppath);
        copyHeaders(hdrs,[cdir '/Solvers/Source/Include/Ooqp/']);
        %Fix couple problems
        str = fileread([ooqppath '/src/Readers/hash.C']);
        str = regexprep(str,'hash(','hashL(');
        fp = fopen([ooqppath '/src/Readers/hash.C'],'w');
        fprintf(fp,'%s',str); fclose(fp);
        str = fileread([ooqppath '/src/QpSolvers/MehrotraSolver.C']);
        if(isempty(strfind(str,'extern double gmu;')))
            str = regexprep(str,'double gmu;','extern double gmu;');
            fp = fopen([ooqppath '/src/QpSolvers/MehrotraSolver.C'],'w');
            fprintf(fp,'%s',str); fclose(fp);
        end
        n=1;
        sdir = [ooqppath '\src'];
        name = 'libooqp';
        opts = [];  
        opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};        
        opts.exclude = {'QpGenSparseOblio.C','QpGenSparseSuperLu.C','OoqpPetscMonitor.C',...
                        'QpGenSparseGondzioDriver.C','QpGenSparseMa57GondzioDriver.C','QpGenSparseMehrotraDriver.C',...
                        'QpBoundDenseGondzioDriver.C','QpGenDenseGondzioDriver.C','SvmGondzioDriver.C','HuberGondzioDriver.C'};
        opts.exFilter = {'QpBoundPetsc*'};
        opts.exFolder = {'Ampl','CInterface','PetscLinearAlgebra','QpExample','Mex'};
        opts.compileAsCpp = true;
        opts.ExcepwCExtern = true;
        %Copy over OPTI modifications
        unzip([cd '/Solvers/Source/ooqp/Ma57Solver.zip'],[cd '/Solvers/Source/ooqp/Ma57']);
        unzip([cd '/Solvers/Source/ooqp/PardisoSolver.zip'],[cd '/Solvers/Source/ooqp/Pardiso']);
        pause(0.1); rehash;
        if(~isempty(MA57))
            copyfile([cd '/Solvers/Source/ooqp/Ma57/Ma57Solver.c'],[ooqppath '\src\LinearSolvers\Ma57Solver\Ma57Solver.c'],'f');
            copyfile([cd '/Solvers/Source/ooqp/Ma57/Ma57Solver.h'],[ooqppath '\src\LinearSolvers\Ma57Solver\Ma57Solver.h'],'f');
            if(strcmpi(MA57,'MATLAB'))
                opts.ex64PP = 'MA57_64BIT_INT'; %64bit integers in API, OPTI MA57 compiled with 32bit ints
                %Check libmwma57 exists, build if required
                checkLibMWMA57();
            end
        end
        if(~isempty(PARDISO))
            if(~exist([ooqppath '\src\LinearSolvers\Pardiso'],'dir')), mkdir([ooqppath '\src\LinearSolvers\'],'Pardiso'); end
            copyfile([cd '/Solvers/Source/ooqp/Pardiso/PardisoSolver.cpp'],[ooqppath '\src\LinearSolvers\Pardiso\PardisoSolver.cpp'],'f');
            copyfile([cd '/Solvers/Source/ooqp/Pardiso/PardisoSolver.h'],[ooqppath '\src\LinearSolvers\Pardiso\PardisoSolver.h'],'f');
            copyfile([cd '/Solvers/Source/ooqp/Pardiso/QpGenSparsePardiso.cpp'],[ooqppath '\src\QpGen\QpGenSparsePardiso.cpp'],'f');
            copyfile([cd '/Solvers/Source/ooqp/Pardiso/QpGenSparsePardiso.h'],[ooqppath '\src\QpGen\QpGenSparsePardiso.h'],'f');
            if(strcmpi(PARDISO,'mkl'))
                opts.mkl = true; %required for mkl pardiso
            end
        end
        rmdir([cd '/Solvers/Source/ooqp/Pardiso'],'s');
        rmdir([cd '/Solvers/Source/ooqp/Ma57'],'s');
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'libooqp'};  
        comps = {'vc'};
        
    case 'pswarm'
        pswmpath = paths{1};
        %Fix couple problems
        str = fileread([pswmpath '/pswarm.h']);
        if(~isempty(strfind(str,'#define SYS_RANDOM 0')))
            str = regexprep(str,'#define SYS_RANDOM 0','#define SYS_RANDOM 1');
            str = regexprep(str,'int sucpollsteps;','int sucpollsteps;int solveriters;');
            fp = fopen([pswmpath '/pswarm.h'],'w');
            fprintf(fp,'%s',str); fclose(fp);
            %C file changes
            str = fileread([pswmpath '/pswarm.c']);
            str = regexprep(str,'struct Stats stats;','extern struct Stats stats; extern struct Options opt;');
            str = strrep(str,'/* some printf can be done here */','stats.solveriters = iter;');
            str = strrep(str,'struct Options opt = {','struct Options optRUB = { ');
            fp = fopen([pswmpath '/pswarm.c'],'w');
            fprintf(fp,'%s',str); fclose(fp);
        end
        %Copy Header Files
        if(~exist([cdir '/Solvers/Source/Include/Pswarm'],'dir'))
            mkdir([cdir '/Solvers/Source/Include/Pswarm']);
        end
        %PSwarm Headers
        [~,hdrs] = VS_BuildFileList(pswmpath,{'include'});
        copyHeaders(hdrs,[cdir '/Solvers/Source/Include/Pswarm/']);
        n = 1;
        sdir = pswmpath;
        name = 'libpswarm';
        opts = [];
        opts.exPP = {'LINEAR','_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE'};
        opts.exclude = {'showcache.c','user.c','pswarm_main.c','pswarm_py.c','pswarm_r.c','cache.c'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'libpswarm'};  
        comps = {'vc'};
        
    case 'csdp'
        csdppath = paths{1};
        %Copy over OPTI modifications
        unzip([cd '/Solvers/Source/Csdp/CSDP_mex_changes.zip'],[cd '/Solvers/Source/Csdp/CSDPChange']);
        pause(0.1); rehash;
        copyfile([cd '/Solvers/Source/Csdp/CSDPChange/easysdp.c'],[csdppath '\lib\easysdp.c'],'f');
        copyfile([cd '/Solvers/Source/Csdp/CSDPChange/sdp.c'],[csdppath '\lib\sdp.c'],'f');
        copyfile([cd '/Solvers/Source/Csdp/CSDPChange/declarations.h'],[csdppath '\include\declarations.h'],'f');
        rmdir([cd '/Solvers/Source/Csdp/CSDPChange'],'s');
        %Copy Header Files
        if(~exist([cdir '/Solvers/Source/Include/Csdp'],'dir'))
            mkdir([cdir '/Solvers/Source/Include/Csdp']);
        end
        %PSwarm Headers
        [~,hdrs] = VS_BuildFileList([csdppath '\include']);
        copyHeaders(hdrs,[cdir '/Solvers/Source/Include/Csdp/']);
        n = 1;        
        sdir = [csdppath '\lib'];
        hdrs = [csdppath '\include'];
        name = 'libcsdp';
        opts.exPP = {'_CRT_SECURE_NO_WARNINGS','USEOPENMP','NOSHORTS'};
        opts.openMP = true;
        opts.exclude = 'user_exit.c';
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'libcsdp'};  
        comps = {'vc'};
        
%           a) In order to accommodate more features of CSDP from MEX I modified 
%           sdp.c and easysdp.c. Changes include (for interest):
%               - easysdp.c:
%                   - Adding struct paramstruc params as input arg
%                   - Set printlevel = 0;
%                   - Adding ppinf,pdinf,prealgap,pxzgap as extra return args
%                   - Modified code to accept, return above appropriately
%               - sdp.c:
%                   - Changed user_exit call to accept any non zero return arg
%                   - Replaced all exit() calls with return instead
%                   - Added two extra return codes for C and constraint checks
%               - declarations.h
%                   - Updated function prototypes as per the above
%                   - Declared malloc to _aligned_malloc so we can align memory for
%                   SSE2 and BLAS/LAPACK (16 byte alignment)
        
    case 'scip'
        error('Not implemented');
        
    case 'asl'        
        aslpath = paths{1};
        %Make required changes
        if(exist([aslpath '/arith.h0'],'file'))
            copyfile([aslpath '/arith.h0'],[aslpath '/arith.h'],'f');
        end
        if(exist([aslpath '/stdio1.h0'],'file'))
            copyfile([aslpath '/stdio1.h0'],[aslpath '/stdio1.h'],'f');
        end
        str = fileread([aslpath '/mainexit.c']);
        if(isempty(strfind(str,'//exit(n);')))
            str = strrep(str,'exit(n);','//exit(n);');
            fp = fopen([aslpath '/mainexit.c'],'w'); fprintf(fp,'%s',str); fclose(fp);
            str = fileread([aslpath '/stderr.c']);
            str = strrep(str,'AllocConsole();','//AllocConsole();');
            fp = fopen([aslpath '/stderr.c'],'w'); fprintf(fp,'%s',str); fclose(fp);
            str = fileread([aslpath '/jac0dim.c']);
            str = strrep(str,sprintf('what_prog();\n\t\tfprintf(Stderr,\n'),'//what_prog();');
            str = regexprep(str,'"jacdim: got M = %d,','//"jacdim: got M = %d,');
            str = strrep(str,'exit(1);','//exit(1);');
            str = strrep(str,'if (n_con < 0 || n_var <= 0 || n_obj < 0) {','if (n_con < 0 || n_var <= 0 || n_obj < 0) { return NULL;');
            fp = fopen([aslpath '/jac0dim.c'],'w'); fprintf(fp,'%s',str); fclose(fp);
        end
        %Copy header files
        if(~exist([cdir '/Utilities/Source/Include/Asl'],'dir'))
            mkdir([cdir '/Utilities/Source/Include/Asl']);
        end
        %Asl Headers
        [~,hdrs] = VS_BuildFileList(aslpath);
        copyHeaders(hdrs,[cdir '/Utilities/Source/Include/Asl/']);
        copyfile([aslpath '/opcode.hd'],[cdir '/Utilities/Source/Include/Asl/opcode.hd'],'f');
        n=1;
        sdir = aslpath;
        name = 'libasl';
        opts = [];
        opts.exPP = {'Arith_Kind_ASL=1','Sscanf=sscanf','Printf=printf','Sprintf=sprintf',...
                     'Fprintf=fprintf','snprintf=_snprintf','NO_STDIO1','_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE'};
        opts.exclude = {'arithchk.c','atof.c','b_search.c','dtoa.c','fpinit.c','funcadd.c',...
                        'funcadd0.c','funcaddk.c','funcaddr.c','obj_adj0.c','sjac0dim.c',...
                        'sprintf.c','sscanf.c','printf.c','mpec_adj0.c'};
        opts.charset = 'MultiByte';            
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'libasl'};  
        comps = {'vc'};
        
    case 'rmathlib'
        error('not supported');

    case 'ma27'
        n = 1;
        sdir = paths{1};
        name = 'libma27';
        opts = [];       
        opts.cpp = false; opts.ifortver = ifortver;
        opts.include = {'ma27.f'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'libma27'};  
        comps = {'if'};
        
    case 'ma57'
        if(length(paths) > 1) %assume metis path 2nd
            opti_VSBuild('metis',paths{2}); %build METIS first
        end
        n = 1;
        sdir = [paths{1} '\src'];
        name = 'libma57';
        opts = [];       
        opts.cpp = false; opts.ifortver = ifortver;
        opts.include = {'ma57d.f','ddeps.f'};
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'libma57'};  
        comps = {'if'};
        
    case 'blas'
        n = 1;
        sdir = paths{1};
        name = 'libblas';
        opts = [];       
        opts.cpp = false; opts.ifortver = ifortver;
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'libblas'};  
        comps = {'if'};
        
    case 'lapack'
        n = 1;
        sdir = paths{1};
        name = 'liblapack';
        opts = [];       
        opts.cpp = false; opts.ifortver = ifortver;
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'liblapack'};  
        comps = {'if'};    
        
    case 'mwma57'
        paths = [cd '\Solvers\Source\opti\libmwma57'];
        unzip([cd '/Solvers/Source/opti/libmwma57.zip'],paths);        
        n = 1;
        sdir = paths;
        name = 'libmwma57';
        opts = [];       
        opts.dll = true; opts.def = 'libmwma57.def';
        VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
        %Write the Solution File
        solpath = VS_WriteSol(VSPRJ,vsver);
        %List of projects to compile and move
        projs = {'libmwma57'};  
        comps = {'vc'}; 
        
    otherwise
        error('Unknown solver!');
end

%Compile Projects             
if(~isempty(projs))
    try
        compileProjects(vsdir,projs,comps,solpath,vsver);
        %Copy out compiled libraries
        copyLibs(solpath,projs,comps,cdir);
        %Special case of mwma57, delete unzipped dir
        if(strcmpi(projs{1},'libmwma57'))
            pause(0.1); rehash;
            rmdir([cd '/Solvers/Source/opti/libmwma57'],'s');
        end
    catch ME
        str = 'There was a problem compiling or moving libraries, please perform this step manually:';
        str = sprintf('%s\n - Compile the following libraries in Visual Studio (Win32 + Win64):\n',str);
        for i = 1:length(projs)
            str = sprintf('%s    - %s\n',str,projs{i});
        end
        str = sprintf('%s\nThen move each of the above libraries (Win32 + Win64) to the following folder:\n',str);
        str = sprintf('%s - %s/Source/libs/win32 or win64\n\n',str,strrep(cdir,filesep,'/'));
        str = sprintf('%sThe Visual Studio Solution is located at:\n - %s\n',str,strrep(solpath,filesep,'/'));
        error([str '\n\nError: %s'],ME.message);
    end
end


function copyLibs(solpath,projs,comps,cdir)

fprintf('Copying Compiled Libraries:\n');
idx = strfind(solpath,filesep);
spath = solpath(1:idx(end));

if(any(strcmp(projs{1},{'libasl','librmathlib'})))
    d = 'Utilities';
else
    d = 'Solvers';
end

for i = 1:length(projs)
    switch(comps{i})
        case {'vc','if'}
            w32lib = [spath 'Release' filesep projs{i} '.lib'];
            w64lib = [spath 'x64' filesep 'Release' filesep projs{i} '.lib'];
        case 'ic'
            w32lib = [spath '..\' projs{i} '\Release' filesep projs{i} '.lib'];
            w64lib = [spath '..\' projs{i} '\x64' filesep '\Release' filesep projs{i} '.lib'];
    end
    fprintf('Copying win32 and win64 libraries of ''%s''...',projs{i}); 
    try
        copyfile(w32lib,[cdir '/' d '/Source/lib/win32/' projs{i} '.lib'],'f');
        ok32 =1;
    catch ME
        fprintf(2,'\nCould not find 32-bit library - ensure it compiled without error. [Error: %s]\n',ME.message);
        ok32 = 0;
    end
    try
        copyfile(w64lib,[cdir '/' d '/Source/lib/win64/' projs{i} '.lib'],'f');
        ok64 = 1;
    catch ME
        fprintf(2,'\nCould not find 64-bit library - ensure it compiled without error. [Error: %s]\n',ME.message);
        ok64 = 0;
    end
    if(ok32&&ok64)
        fprintf('Done!\n');
    end
end
fprintf('\n');


function compileProjects(vsdir,projs,comps,solpath,vsver)

fprintf('Compiling Visual Studio Projects using %s [This May Take a Few Minutes]:\n',upper(vsver));
ccdir = cd;
cd(vsdir);
%estr = ['!devenv ' solpath ' /build Release /project ']; 
n = 1;
for i = 1:length(projs)
    switch(comps{i})
        case {'vc','if'}
            if(strcmpi(comps{i},'vc'))
                cmp = 'Visual C++';
            else
                cmp = 'Intel Fortran';
            end
            if(strcmpi(vsver,'vs2013'))
                estr32 = ['!devenv "' solpath '" /build "Release|Win32" /project ' projs{i} ' /projectconfig "Release|Win32"'];
                estr64 = ['!devenv "' solpath '" /build "Release|x64" /project ' projs{i} ' /projectconfig "Release|x64"'];
            else
                estr32 = ['!devenv "' solpath '" /build Release /project ' projs{i} ' /projectconfig Release|Win32'];
                estr64 = ['!devenv "' solpath '" /build Release /project ' projs{i} ' /projectconfig Release|x64'];
            end
        case 'ic'
            cmp = 'Intel C++';
            idx = strfind(solpath,filesep); %assume not main project
            spath = solpath(1:idx(end-1));
            estr32 = ['!devenv "' spath projs{i} filesep projs{i} '.vcxproj" /build Release|Win32'];
            estr64 = ['!devenv "' spath projs{i} filesep projs{i} '.vcxproj" /build Release|x64'];
    end
    fprintf('Compiling Project (%d of %d) ''%s'' [%s: Win32]...',n,2*length(projs),projs{i},cmp); n = n + 1;
    eval(estr32);
    fprintf('Done!\nCompiling Project (%d of %d) ''%s'' [%s: Win64]...',n,2*length(projs),projs{i},cmp); n = n + 1;
    eval(estr64);
    fprintf('Done!\n');
end
cd(ccdir); 
fprintf('\n');


function copyHeaders(hdrs,dpath)

try
    for i = 1:size(hdrs,1)
        for j = 1:length(hdrs{i,2})
            copyfile([hdrs{i,1} '/' hdrs{i,2}{j}],[dpath hdrs{i,2}{j}],'f');
        end
    end
catch
    fprintf(2,'There was an error copying the header files for the selected solver - you will have to manually copy ALL Solver .h/.hpp files to OPTI/Solvers/Source/Include/$SOLVERNAME$\n');
end


function checkLibMWMA57
%Check to see if libmwma57 exists (pretend library file to allow linking against MATLAB's MA57)
p = [cd '/Solvers/Source/lib'];
switch(computer)
    case 'PCWIN'
        p = [p '/win32/libmwma57.lib'];
    case 'PCWIN64'
        p = [p '/win64/libmwma57.lib'];
end
if(~exist(p,'file'))
    opti_VSBuild('mwma57'); %build it
end

function checkLibMUMPS(path)
%Check to see if libmwma57 exists (pretend library file to allow linking against MATLAB's MA57)
p = [cd '/Solvers/Source/lib'];
switch(computer)
    case 'PCWIN'
        p = [p '/win32/libdmumps_c.lib'];
    case 'PCWIN64'
        p = [p '/win64/libdmumps_c.lib'];
end
if(~exist(p,'file'))
    if(length(path) > 1)
        opts.expaths = path{2};
    else
        error('OPTI Tried to build MUMPS automatically for you, but you need to specify the METIS path when calling opti_VSBuild');
    end
    opti_VSBuild('MUMPS',path,opts); %build it
end

% %% OLD CODE
% BONMIN CODE
% if(haveCPLEX)    
%     [~,cplx_inc] = opti_FindCplex();
%     %Bonmin with CPLEX
%     VSPRJ(n).sdir = VSPRJ(n-1).sdir; VSPRJ(n).hdrs = VSPRJ(n-1).hdrs; VSPRJ(n).opts=VSPRJ(n-1).opts; 
%     VSPRJ(n).name = 'libbonmincplex';
%     VSPRJ(n).hdrs = [VSPRJ(n).hdrs cplx_inc];
%     VSPRJ(n).opts.exPP = [VSPRJ(n).opts.exPP 'COIN_HAS_CPX'];
%     n = n + 1;    
%     % OsiCpx (must add Osi folder manually to include)
%     sdir = [bminpath '\..\Osi\src\OsiCpx'];
%     hdrs = {[bminpath '\..\CoinUtils\src'],[bminpath '\..\BuildTools\headers'],cplx_inc};
%     name = 'libosicpx';
%     opts = [];
%     opts.exPP = {'COINUTILS_BUILD','_CRT_SECURE_NO_WARNINGS'};          
%     VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% end
% if(buildLocals)
%     clppath = 'full path to CLP here'; %e.g. 'C:\Solvers\CLP'
%     cbcpath = 'full path to CBC here'; %e.g. 'C:\Solvers\CBC'
%     ippath = 'full path to IPOPT here'; %e.g. 'C:\Solvers\IPOPT'
%     mumpspath = 'full path to MUMPS here'; %e.g. 'C:\Solvers\MUMPS'
%     %CLP
%     sdir = [clppath '\src'];
%     hdrs = {[clppath '\..\CoinUtils\src'], [clppath '\..\Osi\src'], [clppath '\..\BuildTools\headers']};
%     name = 'libclp';
%     opts = [];
%     opts.exPP = {'CLP_BUILD','_CRT_SECURE_NO_WARNINGS'};
%     opts.exclude = {'MyEventHandler.cpp','MyMessageHandler.cpp','unitTest.cpp','CbcOrClpParam.cpp',...
%                     'ClpCholeskyMumps.cpp','ClpCholeskyUfl.cpp','ClpCholeskyWssmp.cpp','ClpCholeskyWssmpKKT.cpp'};
%     opts.exFilter = {'Abc*','CoinAbc*'};
%     VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
%     % CoinUtils
%     sdir = [clppath '\..\CoinUtils\src'];
%     hdrs = {[clppath '\..\BuildTools\headers']};
%     name = 'libcoinutils';
%     opts = [];
%     opts.exPP = {'COINUTILS_BUILD','_CRT_SECURE_NO_WARNINGS'};          
%     VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
%     % Osi
%     sdir = [clppath '\..\Osi\src'];
%     hdrs = {[clppath '\..\CoinUtils\src'] [clppath '\..\BuildTools\headers']};
%     name = 'libosi';
%     opts = [];
%     opts.exPP = {'OSI_BUILD','_CRT_SECURE_NO_WARNINGS'};   
%     opts.exFolder = {'OsiCpx','OsiGlpk','OsiGrb','OsiMsk','OsiSpx','OsiXpr','OsiCommonTest'};
%     VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
%     % CBC
%     sdir = [cbcpath '\src'];
%     hdrs = {[cbcpath '\..\CoinUtils\src'], [cbcpath '\..\Cgl\src'], [cbcpath '\..\Clp\src'], [cbcpath '\..\Osi\src'], [cbcpath '\..\BuildTools\headers']};
%     name = 'libcbc';
%     opts = [];
%     opts.exPP = {'CBC_BUILD','COIN_FAST_CODE','CLP_FAST_CODE','USE_CBCCONFIG','COIN_NO_TEST_DUPLICATE','_CRT_SECURE_NO_WARNINGS'};
%     opts.exclude = {'unitTest.cpp','unitTestClp.cpp','CoinSolve.cpp','CbcGeneric.cpp','CbcBranchBase.cpp'};
%     VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
%     % CGL
%     sdir = [cbcpath '\..\Cgl\src'];
%     hdrs = {[cbcpath '\..\CoinUtils\src'], [cbcpath '\..\Clp\src\'], [cbcpath '\..\Osi\src'], [cbcpath '\..\BuildTools\headers']};
%     name = 'libcgl';
%     opts = [];
%     opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};          
%     VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
%     % IPOPT (remember to copy updated IpMa57TSolverInterface files and define WIN64 for x64 build)
%     sdir = [ippath '\src'];
%     hdrs = {[mumpspath '\include'], [mumpspath '\libseq'], [ippath '\..\BuildTools\headers']};
%     name = 'libipopt';
%     opts = [];
%     opts.exPP = {'IPOPT_BUILD','FUNNY_MA57_FINT','WIN64','COINHSL_HAS_MA57','_CRT_SECURE_NO_WARNINGS'};
%     opts.exFolder = {'Apps','Algorithm\Inexact','contrib\LinearSolverLoader'};
%     opts.exclude = 'AmplTNLP.cpp';
%     VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;    
%     % BONMIN Solver (optional)
%     sdir = [bminpath '\src'];
%     hdrs = {[bminpath '\..\Cgl\src'],[bminpath '\..\Cbc\src'],[bminpath '\..\Clp\src'],[bminpath '\..\CoinUtils\src'],[bminpath '\..\Ipopt\src'],[bminpath '\..\Osi\src'], [bminpath '\..\BuildTools\headers'],[bminpath '\Examples\CppExample']};
%     name = 'Bonmin';
%     opts = [];
%     opts.exPP = {'BONMIN_BUILD','IPOPT_BUILD','_CRT_SECURE_NO_WARNINGS'};
%     opts.include = {'..\Examples\CppExample\MyBonmin.cpp';'..\Examples\CppExample\MyTMINLP.cpp'};
%     opts.console = true;
%     opts.linklib = {'libbonmin.lib','libipopt.lib','libcbc.lib','libclp.lib','libcgl.lib','libcoinutils.lib','libosi.lib',...
%                     'libdmumps_c.lib','libdmumps_f.lib','libseq_c.lib','libseq_f.lib','libmetis.lib','libpord.lib'};
%     opts.linkpath = {'..\libbonmin','C:\Users\Jonathan Currie\Documents\AUT\Code and Models\Matlab\OPTI\Solvers\mumps\Source\lib\win32'};
%     opts.mkllink = true;
%     VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% end

