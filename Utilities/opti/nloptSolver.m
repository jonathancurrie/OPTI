function [e,ineq,eq,uncon,glob,deriv,subprob] = nloptSolver(varargin)
%NLOPTSOLVER  Return enum for a NLOPT Solver
%
%   [e,ineq,eq] = nloptSolver(name) returns the enumeration e, whether
%   the solver allows inequalities (ineq) and equality (eq) nonlinear
%   constraints.
%
%   nloptSolver() prints a list of all enabled NLOPT Solvers and their
%   functionality.

%   Copyright (C) 2011 Jonathan Currie (I2C2)

global solvers

solvers = cell(31,1);

%Enum Name, Name, G/L, N/D, ineq, eq, infbnd, needsub, enum
solvers{1} = {'AUGLAG',                 'Augmented Lagrangian - Top Layer',                            'G','A',1,1,1,1,36}; i = 2;
% solvers{i} = {'AUGLAG_EQ',              'Augmented Lagrangian for Equalities',                         '','',0,0,0,1,37}; i = i+1; %Not implementing equal AL, local inequal solvers atm
% solvers{i} = {'G_MLSL',                 'Multi-Level Single-Linkage, Random',                           'G','',0,0,0,1,38}; i = i+1; %Seems too slow?
% solvers{i} = {'G_MLSL_LDS',             'Multi-Level Single-Linkage, Quasi-Random',                     'G','',0,0,0,1,39}; i = i+1;
% solvers{i} = {'GD_MLSL',                'Multi-Level Single-Linkage, Random',                           'G','D',0,0,0,0,21}; i = i+1; %Doesn't seem to use gradient?
% solvers{i} = {'GD_MLSL_LDS',            'Multi-Level Single-Linkage, Quasi-Random',                     'G','D',0,0,0,0,23}; i = i+1; 
solvers{i} = {'GD_STOGO',               'StoGO',                                                        'G','D',0,0,0,0,8}; i = i+1;
solvers{i} = {'GD_STOGO_RAND',          'StoGO, Randomized Search',                                     'G','D',0,0,0,0,9}; i = i+1;
solvers{i} = {'GN_CRS2_LM',             'Controlled Random Search, Local Mutation',                     'G','N',0,0,0,0,19}; i = i+1;
solvers{i} = {'GN_DIRECT',              'DIviding RECTangles',                                          'G','N',0,0,0,0,0}; i = i+1;
solvers{i} = {'GN_DIRECT_NOSCAL',       'Unscaled DIviding RECTangles',                                 'G','N',0,0,0,0,3}; i = i+1;
solvers{i} = {'GN_DIRECT_L',            'DIviding RECTangles-L, Biased Local',                          'G','N',0,0,0,0,1}; i = i+1;
solvers{i} = {'GN_DIRECT_L_NOSCAL',     'Unscaled DIviding RECTangles-L, Biased Local',                 'G','N',0,0,0,0,4}; i = i+1;
solvers{i} = {'GN_DIRECT_L_RAND',       'Randomized DIviding RECTangles-L, Biased Local',               'G','N',0,0,0,0,2}; i = i+1;
solvers{i} = {'GN_DIRECT_L_RAND_NOSCAL','Unscaled, Randomized DIviding RECTangles-L, Biased Local',     'G','N',0,0,0,0,5}; i = i+1;
solvers{i} = {'GN_ORIG_DIRECT',         'Original DIviding RECTangles, Gablonsky Imp.',                 'G','N',0,0,0,0,6}; i = i+1; %these appear to crash with ineq con
solvers{i} = {'GN_ORIG_DIRECT_L',       'Original DIviding RECTangles-L, Biased Local, Gablonsky Imp.', 'G','N',0,0,0,0,7}; i = i+1;
solvers{i} = {'GN_ISRES',               'Improved Stochastic Ranking Evolution Strategy',               'G','N',1,1,0,0,35}; i = i+1;
solvers{i} = {'GN_ESCH',                'ESCH Evolutionary Algorithm',                                  'G','N',0,0,0,0,42}; i = i+1;
solvers{i} = {'GN_MLSL',                'Multi-Level Single-Linkage, Random',                           'G','N',0,0,0,0,20}; i = i+1;
solvers{i} = {'GN_MLSL_LDS',            'Multi-Level Single-Linkage, Quasi-Random',                     'G','N',0,0,0,0,22}; i = i+1;

solvers{i} = {'LD_AUGLAG',              'Augmented Lagrangian - Top Layer',                            'L','D',1,1,1,1,31}; i = i+1;
% solvers{i} = {'LD_AUGLAG_EQ',           'Augmented Lagrangian for Equalities',                         'L','D',1,1,0,1,33}; i = i+1;
solvers{i} = {'LD_CCSAQ',               'Conservative Convex Separable Approximation Quadratic',        'L','D',0,0,1,0,41}; i = i+1;
solvers{i} = {'LD_LBFGS',               'Limited Memory BFGS',                                          'L','D',0,0,1,0,11}; i = i+1;
% solvers{i} = {'LD_LBFGS_NOCEDAL',       'Limited Memory BFGS (Nocedal Version)',           'L','D',0,0,0,10}; i = i+1; %Not available in current release
solvers{i} = {'LD_MMA',                 'Method of Moving Asymptotes',                                  'L','D',0,0,1,0,24}; i = i+1; %Note nl con didn't seem to work
solvers{i} = {'LD_TNEWTON',             'Truncated Newton',                                             'L','D',0,0,1,0,15}; i = i+1;
solvers{i} = {'LD_TNEWTON_RESTART',     'Truncated Newton with Restarting',                             'L','D',0,0,1,0,16}; i = i+1;
solvers{i} = {'LD_TNEWTON_PRECOND',     'Preconditioned Truncated Newton',                              'L','D',0,0,1,0,17}; i = i+1;
solvers{i} = {'LD_TNEWTON_PRECOND_RESTART','Preconditioned Truncated Newton with Restarting',           'L','D',0,0,1,0,18}; i = i+1;
solvers{i} = {'LD_SLSQP',               'Sequential Least Squares Quadratic Programming',               'L','D',1,1,1,0,40}; i = i+1;
solvers{i} = {'LD_VAR1',                'Limited Memory Variable-Metric, Rank 1',                       'L','D',0,0,1,0,13}; i = i+1;
solvers{i} = {'LD_VAR2',                'Limited Memory Variable-Metric, Rank 2',                       'L','D',0,0,1,0,14}; i = i+1;
solvers{i} = {'LN_AUGLAG',              'Augmented Lagrangian - Top Layer',                            'L','N',1,1,1,1,30}; i = i+1;
% solvers{i} = {'LN_AUGLAG_EQ',           'Augmented Lagrangian for Equalities',                       'L','N',1,1,0,1,32}; i = i+1;
solvers{i} = {'LN_BOBYQA',              'Bound-Constrained Optimization via Quadratic Models',          'L','N',0,0,1,0,34}; i = i+1;
solvers{i} = {'LN_COBYLA',              'Constrained Optimization by Linear Approximations',            'L','N',1,1,1,0,25}; i = i+1;
solvers{i} = {'LN_NELDERMEAD',          'Nelder-Mead Simplex',                                          'L','N',0,0,1,0,28}; i = i+1;
solvers{i} = {'LN_NEWUOA',              'Unconstrained Optimization via Quadratic Models',              'L','N',0,0,1,0,26}; i = i+1;
% solvers{i} = {'LN_NEWUOA_BOUND',        'Bound Constrained Optimization via NEWUOA Quadratic Models',   'L','N',0,0,0,0,27}; i = i+1; %Didnt seem to work
solvers{i} = {'LN_PRAXIS',              'Principle Axis',                                               'L','N',0,0,1,0,12}; i = i+1;
solvers{i} = {'LN_SBPLX',               'Subplex variation of Nelder-Mead',                             'L','N',0,0,1,0,29}; i = i+1;

if(nargin==0 && nargout ==0)
    printSolvers();
    return;
end

arg = varargin{1};

if(ischar(arg))
    j = findByName(arg);
    e = solvers{j}{9}; %enum
elseif(isnumeric(arg))
    j = findByEnum(arg);
    e = solvers{j}{1}; %name
else
    error('You can only search for a solver by name or enum');
end
%Assign common props (ineq,eq,uncon,glob,deriv)
ineq = solvers{j}{5};
eq = solvers{j}{6};
uncon = solvers{j}{7};
glob = strcmp(solvers{j}{3},'G');
deriv = strcmp(solvers{j}{4},'D');
subprob = solvers{j}{8};


%Search Helper Functions
function j = findByName(arg)
global solvers

for i = 1:length(solvers)
    if(strcmpi(arg,solvers{i}{1}))
        j = i;
        return;
    end
end
%If here, we've failed
error('Unknown NLOPT solver: %s',arg);

function j = findByEnum(arg)
global solvers

for i = 1:length(solvers)
    if(arg == solvers{i}{9})
        j = i;
        return;
    end
end
%If here, we've failed
error('Unknown NLOPT solver enum: %s',arg);


function printSolvers()

global solvers

len = length(solvers);

fprintf('\n-----------------------------------------------------------------------------------------------------------------------------------------------\n')
fprintf('NLOPT AVAILABLE SOLVERS:\n\n');

fprintf(' NAME                          SCOPE     DERIV  UNCON  INEQ    EQ    DESCRIPTION                                                    SUBOPT\n');  
for i = 1:len
    fprintf('%-28s | %-8s | %-4s | %-4s | %-4s | %-4s | %-60s | %-5s |\n',solvers{i}{1},getScope(solvers{i}{3}),getReq(solvers{i}{4}),getYes(solvers{i}{7}),...
            getYes(solvers{i}{5}),getYes(solvers{i}{6}),solvers{i}{2},getReqD(solvers{i}{8}));
end

fprintf('\n-----------------------------------------------------------------------------------------------------------------------------------------------\n')

%Print Helper Functions
function str = getScope(g)
switch(g)
    case 'G'
        str = 'Global';
    case 'L'
        str = 'Local';
    case 'A'
        str = 'Any';
end

function str = getReq(g)
switch(g)
    case {'D'}
        str = 'Req';
    case {'N'}
        str = '';
    case 'A'
        str = 'Any';
end

function str = getReqD(g)
if(g)
    str = 'Req';
else
	str = '';
end

function str = getYes(u)
if(u)
    str = 'YES';
else
    str = '';
end


