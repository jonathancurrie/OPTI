function [x,info] = opti_mumps(A,b)
%OPTI_MUMPS Solve a system of sparse linear equations using MUMPS
%
%   x = opti_mumps(A,b) solves the system of sparse linear equations Ax = b.
%
%   THIS IS A WRAPPER FOR MUMPS USING THE MEX INTERFACE
%   See supplied License

%   Copyright (C) 2011 Jonathan Currie (IPL)

if(~isreal(A) || ~isreal(b))
    error('opti_mumps only solves real problems. Use opti_zmumps for complex problems.');
end

t = tic;

%Initialize MUMPS
id = struct('SYM',0,'JOB',-1,'ICNTL',zeros(1,40)-9998,'CNTL',zeros(1,15)-9998,...
            'PERM_IN',-9999,'COLSCA',-9999,'ROWSCA',-9999,'RHS',-9999,'INFOG',zeros(1,40)-9998,...
            'RINFOG',zeros(1,40)-9998,'VAR_SCHUR',-9999,'SCHUR',-9999,'INST',-9999,'SOL',-9999,...
            'REDRHS',-9999,'PIVNUL_LIST',-9999,'MAPPING',-9999,'SYM_PERM',-9999,'UNS_PERM',-9999,'TYPE',0);

id = dmumps(id);

id.Type = 1; %arithmetic type
id.JOB = 6; %analyze + factorize + solve

%Setup and Solve
id.RHS = b;
if(~issparse(A))
    optiwarn('opti:sparse','The A matrix should be sparse, correcting: [sparse(A)]');
    A = sparse(A);
end
id = dmumps(id,A);
x = id.SOL;

%Destroy mumps instance
id.JOB = -2;
id = dmumps(id); %#ok<NASGU>

info.Iterations = [];
info.Time = toc(t);
info.Algorithm = 'MUMPS: Sequential Build';
info.Status = [];
