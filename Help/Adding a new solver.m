%% Steps for Adding a New Solver

% 1) Open checkSolver.m add to 
%       - problem type list
%       - description list
%       - known solver switch statement
% 2) Open optiSolverInfo.m add an entry for the solver and add:
%       - config information
%       - constraint information 
%       - derivative information
%       - unique details
% 3) Open buildConfig.m add to
%       - known solver switch statement
%       - build configuration function
% 4) Open solveOpti.m add to
%       - problem types which can be solved
% 5) Open opti_MEX_Install add install file
% 6) Open opti_Dist_Test add testing file