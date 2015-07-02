%% SCIP Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile Solving Constraint Integer Programs (SCIP) 
% for use with MATLAB. 

% My build platform:
% - Windows 7 x64
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

% %% Visual Studio Builder Commands [LATEST SCIP REQ C++11]
% %NOTE  - the PaxHeader folders change - update as required
% clear
% path = 'C:\Solvers\scipoptsuite-3.2.0\scip-3.2.0'; %e.g. 'C:\Solvers\SCIP'
% splxpath = 'C:\Solvers\scipoptsuite-3.2.0\soplex-2.2.0'; %e.g. 'C:\Solvers\SOPLEX'
% ipoptpath = 'C:\Solvers\Ipopt-3.12.3\Ipopt'; %e.g. 'C:\Solvers\IPOPT'
% n = 1;
% % SCIP
% sdir = [path '\src'];
% cppadpath = [sdir '\cppad']; %use built in version
% hdrs = {[ipoptpath '\src'], [cd '\Solvers\Source\Include\BuildTools'], [splxpath '\src'], cppadpath};
% name = 'libscip';
% opts = [];
% opts.exPP = {'IPOPT_BUILD','_CRT_SECURE_NO_WARNINGS','NO_RAND_R','NO_SIGACTION','NO_STRERROR_R',...
%              'NO_STRTOK_R','ROUNDING_MS','NPARASCIP','WITH_SCIPDEF'};          
% opts.exclude = {'exprinterpret_none.c','nlpi_ipopt_dummy.c','nlpi_xyz.c','lpi_none.c',...
%                 'lpi_clp.cpp','lpi_cpx.c','lpi_grb.c','lpi_msk.c','lpi_qso.c',...
%                 'lpi_spx.cpp','lpi_xprs.c','sorttpl.c','cmain.c',... %check whether we want lpi_spx2 or spx...?
%                 'cppmain.cpp','disp_xyz.c','branch_xyz.c','event_xyz.c','cons_xyz.c',...
%                 'heur_xyz.c','presol_xyz.c','prop_xyz.c','pricer_xyz.c','nodesel_xyz.c',...
%                 'relax_xyz.c','reader_xyz.c','sepa_xyz.c','dialog_xyz.c','compr_xyz.c'};           
% VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% % SOPLEX
% sdir = [splxpath '\src'];
% name = 'libsoplex';
% opts = [];
% opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};
% opts.exclude = {'soplexmain.cpp','simpleexample.cpp'};
% VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% % SCIP executable (to use this, copy all the .libs listed below into the respective .lib directories in libscip)
% sdir = [path '\src'];
% hdrs = {[ipoptpath '\src'], [cd '\Solvers\Source\Include\BuildTools'], [splxpath '\src'], cppadpath};
% name = 'scip';
% opts = [];
% opts.exPP = {'IPOPT_BUILD','_CRT_SECURE_NO_WARNINGS','NO_RAND_R','NO_SIGACTION','NO_STRERROR_R',...
%              'NO_STRTOK_R','ROUNDING_MS','NPARASCIP'};
% opts.include = {'cppmain.cpp'};
% opts.console = true;
% opts.linklib = {'libscip.lib','libsoplex.lib','libipoptbm.lib','libdmumps_c.lib','libdmumps_f.lib','libseq_c.lib','libseq_f.lib','libmetis.lib','libpord.lib','libma57.lib'};
% opts.linkpath = {'..\libscip'};
% opts.mkllink = true;
% VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% %Write the Solution File
% VS_WriteSol(VSPRJ)
% %%

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
%       - FIND #include "lpi/           REPLACE #include "../lpi/
%       - FIND #include "objscip/       REPLACE #include "../objscip/
%       - FIND #include "xml/           REPLACE #include "../xml/
%       - FIND #include "dijkstra/      REPLACE #include "../dijkstra/
%       - FIND #include <cppad/         REPLACE #include <../cppad/
%       - FIND # include <cppad/        REPLACE # include <../cppad/
%   b) In src/cppad/configure.hpp for 32bit change lines 120 and 134
%           # define CPPAD_SIZE_T_SAME_UNSIGNED_INT 1
%           # define CPPAD_TAPE_ADDR_TYPE size_t
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
%       - spxdefines.h from SOPLEX
%   to 
%
%   OPTI/Solvers/Source/Include/Scip

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
 
%MEX Interface Source Files 
src = {'scip/scipmex.cpp scip/scipeventmex.cpp scip/scipnlmex.cpp'};
%Include Directories
inc = {'scip/Include','Include/Scip','Include/scip/nlpi','Include/scip/blockmemshell',...
       'Include/Ipopt'};
%Lib Names [static libraries to link against]
libs = {'libscip','libsoplex','libipoptbm'};
%Options (note options from above used here too)
opts = [];
opts.verb = false;
opts.blas = 'MKL'; 
opts.ma57 = 'Matlab';
opts.mumps = true;
opts.asl = true; %link ASL library

%Compile
opti_solverMex('scip',src,inc,libs,opts);