%% SCIP-SDP Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (IPL)

% This file will help you compile Solving Constraint Integer Programs (SCIP) 
% Semidefinite Program (SDP) Extension for use with MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 0) Complete Compilation as per OPTI instructions for MUMPS and IPOPT, in 
% that order.

% 1) Get SCIP & SoPlex
% SCIP is available from http://scip.zib.de/. Ensure you download the SCIP
% optimization suite, which includes SoPlex. We will create the VS projects 
% below.

% 2) Compile SCIP & SoPlex
% The easiest way to compile SCIP is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:
% 
%% Visual Studio Builder Commands
clear
path = 'C:\Users\Jonathan Currie\Documents\AUT\Code and Models\Optimization Solvers\scipoptsuite-3.0.2\scipoptsuite-3.0.2\scip-3.0.2';%'full path to SCIP here'; %e.g. 'C:\Solvers\SCIP'
sdppath = 'C:\Users\Jonathan Currie\Documents\AUT\Code and Models\Optimization Solvers\scipoptsuite-3.0.2\scipoptsuite-3.0.2\SCIP_SDP';
dsdppath = 'C:\Users\Jonathan Currie\Documents\AUT\Code and Models\Optimization Solvers\DSDP5.8\DSDP5.8';
splxpath = 'full path to SOPLEX here'; %e.g. 'C:\Solvers\SOPLEX'
% ipoptpath = 'full path to IPOPT here'; %e.g. 'C:\Solvers\IPOPT'
n = 1;
% SCIP
sdir = [path '\src'];
cppadpath = [sdir '\cppad']; %use built in version
hdrs = {cppadpath};
name = 'libscip';
opts = [];
opts.exPP = {'IPOPT_BUILD','_CRT_SECURE_NO_WARNINGS','NO_RAND_R','NO_SIGACTION','NO_STRERROR_R',...
             'NO_STRTOK_R','NO_NEXTAFTER','ROUNDING_MS','NPARASCIP'};          
opts.exclude = {'exprinterpret_cppad.cpp','nlpi_ipopt.cpp','nlpi_xyz.c',...
                'lpi_clp.cpp','lpi_cpx.c','lpi_grb.c','lpi_msk.c','lpi_qso.c',...
                'lpi_spx.cpp','lpi_spx121.cpp','lpi_spx132.cpp','lpi_xprs.c','sorttpl.c','cmain.c',...
                'cppmain.cpp','disp_xyz.c','branch_xyz.c','event_xyz.c','cons_xyz.c',...
                'heur_xyz.c','presol_xyz.c','prop_xyz.c','pricer_xyz.c','nodesel_xyz.c',...
                'relax_xyz.c','reader_xyz.c','sepa_xyz.c','dialog_xyz.c'};
opts.exFolder = {'PaxHeaders.29102'};            
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% SOPLEX
sdir = [splxpath '\src'];
name = 'libsoplex';
opts = [];
opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};
opts.exclude = {'soplexmain.cpp','simpleexample.cpp'};
opts.exFolder = {'PaxHeaders.29102'}; 
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% SCIP SDP
sdir = [sdppath '\src'];
hdrs = {[path '\src'], [dsdppath '\include']};
name = 'libscipsdp';
opts = [];
opts.exPP = {'_CRT_SECURE_NO_WARNINGS','USE_DSDP'};
opts.exclude = {'main.cpp'};
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% SCIP-SDP executable
sdir = [sdppath '\src'];
hdrs = {[path '\src']};
name = 'scipsdp';
opts = [];
opts.exPP = {'_CRT_SECURE_NO_WARNINGS','NO_RAND_R','NO_SIGACTION','NO_STRERROR_R',...
             'NO_STRTOK_R','NO_NEXTAFTER','ROUNDING_MS','NPARASCIP'};     
opts.include = {'main.cpp'};
opts.console = true;
opts.mkllink = true;
opts.linklib = {'libscip.lib','libscipsdp.lib','libdsdp.lib'};
opts.linkpath = {'..\..\scip-3.0.2\libscip','C:\Users\Jonathan Currie\Documents\AUT\Code and Models\Matlab\OPTI\Solvers\dsdp\Source\lib\win32'};
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
%Write the Solution File
VS_WriteSol(VSPRJ)
%%

% SCIP SDP changes
% remove scip from include path
% remove objscip from include path
% remove blockmemshell from include path
% remove private declaration line 93 BlockMemoryAllocator.h


% Once complete, you will have a directory called SCIP/libscip. Open the
% Visual Studio 2012 solution file, then complete the following steps:
%   a) (SCIP) Due to the way we have created the project the header files 
%   cannot be found. Perform a find and replace on the CURRENT PROJECT of the
%   following terms:
%       - FIND #include "tclique/       REPLACE #include "../tclique/
%       - FIND #include "blockmemshell/ REPLACE #include "../blockmemshell/
%       - FIND #include <blockmemshell/ REPLACE #include <../blockmemshell/
%       - FIND #include "scip/          REPLACE #include "../scip/
%       - FIND #include "nlpi/          REPLACE #include "../nlpi/
%       - FIND #include "objscip/       REPLACE #include "../objscip/
%       - FIND #include "xml/           REPLACE #include "../xml/
%       - FIND #include "dijkstra/      REPLACE #include "../dijkstra/
%       - FIND #include <cppad/         REPLACE #include <../cppad/
%       - FIND # include <cppad/        REPLACE # include <../cppad/
%   b) In src/cppad/configure.hpp for 32bit change lines 87 and 101 to
%           # define CPPAD_SIZE_T_SAME_UNSIGNED_INT 1
%           # define CPPAD_TAPE_ADDR_TYPE unsigned int
%      while for 64bit change it to
%           # define CPPAD_SIZE_T_SAME_UNSIGNED_INT 0
%           # define CPPAD_TAPE_ADDR_TYPE size_t
%   c) Build a Win32 or x64 Release of each project to compile the code.
%   d) Copy the generated .lib files to the following folder:
%
%   OPTI/Solvers/scip/Source/lib/win32 or win64
%
%   Next we need to copy all the required header files. Copy the header files
%   (and folders) of the following folders:
%       - blockmemshell
%       - nlpi
%       - objscip
%       - scip
%       - cppad (just configure.hpp)
%   to 
%
%   OPTI/Solvers/scip/Source/Include

% 3) MEX Interface
% The MEX interface supplied with SCIP (at the time of this development)
% was quite basic, thus has been updated for use with OPTI. Therefore I
% suggest you use the version supplied with OPTI.

% 4) ASL Interface
% If you wish to be able to solve AMPL .nl models using this interface, you
% can enable this functionality below with "haveASL". However you will need
% to complete building the ASL library, as detailed in Utilities/File
% IO/opti_AMPL_Install.m. Also ensure reader_nl.c and .h are placed in
% scip/Source/AMPL/. Note a bug exists in v3.0.1, so ensure the latest
% release is used.

% 5) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the SCIP MEX file. Once you have completed all the
% above steps, simply run this file to compile SCIP! You MUST BE in the 
% base directory of OPTI!

clear scipsdp

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();
% Dependency Paths
dsdpdir = '..\..\dsdp\Source\';
ipoptdir = '..\..\ipopt\Source\';
mumpsdir = '..\..\mumps\Source\';

fprintf('\n------------------------------------------------\n');
fprintf('SCIPSDP MEX FILE INSTALL\n\n');
   
%Get IPOPT Libraries
post = [' -IInclude -I..\..\ipopt\Source\Include\Common\ -L' ipoptdir libdir ' -llibIpopt']; 
%Optionally Add MA57
post = [post ' -llibmwma57 -DhaveMA57'];
%Get MUMPS Libraries
post = [post ' -I' mumpsdir 'Include -L' mumpsdir libdir '-llibdmumps_c -llibdmumps_f -llibseq_c -llibseq_f -llibmetis -llibpord'];
%Get DSDP Libraries
post = [post ' -I' dsdpdir 'Include -L' dsdpdir libdir '-llibdsdp'];
%Get Intel Fortran Libraries (for MUMPS build) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link];
%Get SCIP Includes and Libraries
post = [post ' -IInclude -IInclude/scipsdp -IInclude/objscip -IInclude/scip -L' libdir ' -llibscipsdp -llibscip -llibsoplex -llibut -output scipsdp'];
   
%CD to Source Directory
cdir = cd;
cd 'Solvers/scip/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims scipsdpmex.cpp scipeventmex.cpp';
try
    eval([pre post])
    movefile(['scipsdp.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:scipsdp','Error Compiling SCIPSDP!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');
