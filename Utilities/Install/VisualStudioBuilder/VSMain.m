%% Visual Studio Project Builder

% Often I found myself manually creating Visual Studio 2010 projects for
% solvers with lots of files, even more directories to include, and having
% to repeat it for each release! This code, together with associated files
% aims to remove the need for this by automatically generating VS2010
% projects!
%
% The code is an afternoon's work and is to be honest, rubbish. But until I
% have a need for 'efficient Visual Studio project creation' it will
% suffice! Minimal error checking so watch out.
%
% Note the code does not help you compile the project, it only assembles
% all the files and includes for you. You will still need to manually set
% the required compilation options.
%
% An example for creating a VS2010 project for BONMIN is listed below.
%
% J.Currie June 2012


%% Example for creating a Visual Studio project for BONMIN
clc
%BONMIN base directory (replace with your solver here)
bdir = 'C:\Users\Jonathan Currie\Documents\AUT\MINLP\bonmin_stable\';

%Desired Project Name
projName = 'libBonmin';

%Bonmin Source Folder (replace with your solver source folder)
src = [bdir 'Bonmin\src'];

%Extra header files to include (optional)
hdrs = {[bdir 'Cgl\src'],[bdir 'Cbc\src'],[bdir 'Clp\src'],[bdir 'CoinUtils\src'],[bdir 'Ipopt\src'],[bdir 'Osi\src']};

%Call the top function to create the project
VS_WriteProj(src,projName,hdrs)

%Open the resulting libBonmin VC++ Project and VS will generate a solution
%for you, easy!

%% Example for SCIP
clc

%Project Name
projName = 'libscip';

%SCIP Source Directory
sdir = 'C:\Users\Jonathan Currie\Documents\AUT\MINLP\scipoptsuite-3.0.1\scip-3.0.1\src';

%Extra Header files to include
hdrs = {'C:\Users\Jonathan Currie\Documents\AUT\MINLP\Ipopt-3.10.2\Ipopt-3.10.2\Ipopt\src'};

%Create the project
VS_WriteProj(sdir,projName,hdrs)

%% Example for DSDP
clc

%Source Directory
sdir = 'C:\Users\Jonathan\Documents\AUT\MINLP\DSDP5.8\DSDP5.8\src';

%Extra Header files to include
hdrs = {'C:\Users\Jonathan\Documents\AUT\MINLP\DSDP5.8\DSDP5.8\include'};

%Project Name
projName = 'libdsdp';

%Create the project
VS_WriteProj(sdir,projName,hdrs)

%% CppAD
clc

%Source Directory
sdir = 'C:\Users\jocurrie\Documents\AUT\Code\Optimization Solvers\cppad-20131125\cppad';

%Project Name
projName = 'libcppad';

%Create the project
VS_WriteProj(sdir,projName,[],[])
