function options = cbcset(varargin)
%CBCSET  Create or alter the options for Optimization with CBC
%
% options = cbcset('param1',value1,'param2',value2,...) creates an CBC
% options structure with the parameters 'param' set to their corresponding
% values in 'value'. Parameters not specified will be set to the CBC
% default.
%
% options = cbcset(oldopts,'param1',value1,...) creates a copy of the old
% options 'oldopts' and then fills in (or writes over) the parameters
% specified by 'param' and 'value'.
%
% options = cbcset() creates an options structure with all fields set to
% cbcset defaults.
%
% cbcset() prints a list of all possible fields and their function.
%
% See supplied CBC Documentation for further details of these options.

%   Copyright (C) 2013 Jonathan Currie (Control Engineering)

% Print out possible values of properties.
if((nargin == 0) && (nargout == 0))
    printfields();
    return
end

Names = {'dualTol','primalTol','allowableGap','allowableFracGap','startalg','strategy','preprocess','heuristics',...
         'nodeStrategy','increment','maxSolutions','strongBranching','trustPseudoCosts','combineSolution','dins',...
         'divingSome','divingCoefficient','divingFractional','divingGuided','divingLineSearch','divingPseudoCost',...
         'divingVectorLength','feasibilityPump','feasibilityPumpPasses','greedyHeuristic','localTreeSearch',...
         'pivotAndFix','randomizedRounding','rens','rins','roundingHeuristic','vubHeuristic','costStrategy',...
         'useCuts','cutDepth','maxCutPassesRoot','maxCutPasses','slowCutPasses','cliqueCuts','flowCoverCuts',...
         'GMICuts','lagomoryCuts','gomoryCuts','knapsackCuts','liftAndProjectCuts','mixedIntegerRoundingCuts',...
         'probingCuts','reduceAndSplitCuts','residualCapacityCuts','twoMirCuts','latwoMirCuts','zeroHalfCuts',...
         'idiotCrash','sprintCrash','crash','factorization','maxFactor','crossover','dualPivot','primalPivot',...
         'dualBound','primalWeight','perturbation','scaling','presolve','preTolerance','passPresolve'};
Defaults = {1e-7,1e-7,0,0,[],1,'SOS','On','Fewest',1e-5,-1,5,5,'On','Off','Off','On','Off','Off','Off','Off',...
            'Off','On',30,'Off','Off','Off','Off','Off','On','On',-1,'Off','On',-1,-1,1,10,'IfMove','IfMove',...
            'Off','Off','IfMove','IfMove','Off','IfMove','On','Off','Off','Root','Off','Off',-1,-1,'Off','Normal',...
            200,'On','automatic','automatic',1e10,1e10,'On','automatic','On',1e-8,5}; 

%Enter and check user args
try
    options = opticheckset(Names,Defaults,@checkfield,varargin{:});
catch ME
    throw(ME);
end


function checkfield(field,value)
%Check a field contains correct data type
switch lower(field)
    %Scalar double
    case {'increment'}
        err = opticheckval.checkScalarDbl(value,field);    
    %Scalar non zero double
    case {'primaltol','dualtol','dualbound','primalweight','pretolerance'}
        err = opticheckval.checkScalarGrtZ(value,field);    
    %Scalar non negative double
    case {'allowablegap','allowablefracgap'}
        err = opticheckval.checkScalarNonNeg(value,field); 
    %Scalar non negative integer
    case {'strongbranching'}
        err = opticheckval.checkScalarIntNonNeg(value,field);
    %Scalar integer -1 or > 0
    case 'maxsolutions'
        err = opticheckval.checkScalarIntM1GrtZ(value,field);
    %Scalar integer with bounds
    case 'strategy'
        err = opticheckval.checkScalarIntBoundLELE(value,field,0,2); 
    case 'trustpseudocosts'
        err = opticheckval.checkScalarIntBoundLELE(value,field,-3,2000000000); 
    case 'feasibilitypumppasses'
        err = opticheckval.checkScalarIntBoundLELE(value,field,0,1e4); 
    case 'vubheuristic'
        err = opticheckval.checkScalarIntBoundLELE(value,field,-2,20); 
    case 'cutdepth'
        err = opticheckval.checkScalarIntBoundLELE(value,field,-1,999999); 
    case {'maxcutpassesroot','maxcutpasses'}
        err = opticheckval.checkScalarIntBoundLELE(value,field,-9999999,9999999); 
    case {'slowcutpasses','idiotcrash','sprintcrash'}
        err = opticheckval.checkScalarIntBoundLEL(value,field,-1,Inf);
    case 'maxfactor'
        err = opticheckval.checkScalarIntBoundLEL(value,field,1,999999);
    case 'passpresolve'
        err = opticheckval.checkScalarIntBoundLEL(value,field,-200,100);
    %String On or Off
    case {'heuristics','combinesolution','dins','localtreesearch','perturbation',}
        err = opticheckval.checkOnOff(value,field);
    %Misc String methods
    case 'startalg'
        err = opticheckval.checkValidString(value,field,{'Primal','Dual','Barrier'});
    case 'preprocess'
        err = opticheckval.checkValidString(value,field,{'Off','On','Equal','SOS','TrySOS','EqualAll','Strategy','Aggregate','ForceSOS'});
    case 'nodestrategy'
        err = opticheckval.checkValidString(value,field,{'Hybrid','Fewest','Depth','UpFewest','DownFewest','UpDepth','DownDepth'});
    case {'divingsome','divingcoefficient','divingfractional','divingguided','divinglinesearch','divingpseudocost','divingvectorlength',...
          'feasibilitypump','greedyheuristic','pivotandfix','randomizedrounding','rins','roundingheuristic'}
        err = opticheckval.checkValidString(value,field,{'Off','On','Both','Before'});
    case 'rens'
        err = opticheckval.checkValidString(value,field,{'Off','On','Both','Before','200','1000','10000','dj','djbefore'});
    case 'coststrategy'
        err = opticheckval.checkValidString(value,field,{'Off','Priorities','ColumnOrder','01first','01last','Length','Singletons','Nonzero'});
    case {'usecuts','liftandprojectcuts','reduceandsplitcuts','residualcapacitycuts'}
        err = opticheckval.checkValidString(value,field,{'Off','On','Root','IfMove','ForceOn'});
    case {'cliquecuts','flowcovercuts','mixedintegerroundingcuts','zerohalfcuts'}
        err = opticheckval.checkValidString(value,field,{'Off','On','Root','IfMove','ForceOn','OnGlobal'});
    case 'knapsackcuts'
        err = opticheckval.checkValidString(value,field,{'Off','On','Root','IfMove','ForceOn','OnGlobal','ForceAndGlobal'});
    case 'gmicuts'
        err = opticheckval.checkValidString(value,field,{'Off','On', 'Root', 'IfMove', 'ForceOn', 'EndOnly', 'Long', 'LongRoot', 'LongIfMove', 'ForceLongOn', 'LongEndOnly'});
    case {'lagomorycuts','latwomircuts'}
        err = opticheckval.checkValidString(value,field,{'Off','EndOnlyRoot', 'EndCleanRoots', 'EndBothRoot', 'EndOnly', 'EndClean', 'EndBoth', 'OnlyAsWell', 'CleanAsWell',...
                                                         'BothAsWell', 'OnlyInstead', 'CleanInstead', 'BothInstead'});
    case 'gomorycuts'
        err = opticheckval.checkValidString(value,field,{'Off', 'On', 'Root', 'IfMove', 'ForceOn', 'OnGlobal', 'ForceAndGlobal', 'ForceLongOn', 'Long'});
    case 'probingcuts'
        err = opticheckval.checkValidString(value,field,{'Off', 'On', 'Root', 'IfMove', 'ForceOn', 'OnGlobal', 'ForceOnGlobal', 'ForceOnBut', 'ForceOnStrong', 'ForceOnButStrong', 'StrongRoot'});
    case 'twomircuts'
        err = opticheckval.checkValidString(value,field,{'Off', 'On', 'Root', 'IfMove', 'ForceOn', 'OnGlobal', 'ForceAndGlobal', 'ForceLongOn'});
    case 'crash'
        err = opticheckval.checkValidString(value,field,{'Off', 'On', 'solow_halim', 'lots'});
    case 'factorization'
        err = opticheckval.checkValidString(value,field,{'Normal','Dense','Simple','OSL'});
    case 'crossover'
        err = opticheckval.checkValidString(value,field,{'On','Off','Maybe','Presolve'});
    case 'dualpivot'
        err = opticheckval.checkValidString(value,field,{'automatic','dantzig','partial','steepest'});
    case 'primalpivot'
        err = opticheckval.checkValidString(value,field,{'automatic', 'exact', 'dantzig', 'partial', 'steepest', 'change', 'sprint'});
    case 'scaling'
        err = opticheckval.checkValidString(value,field,{'automatic', 'equilibrium', 'geometric', 'dynamic', 'rowsonly'});
    case 'presolve'
        err = opticheckval.checkValidString(value,field,{'On','Off','More'});
    otherwise  
        err = MException('OPTI:SetFieldError','Unrecognized parameter name ''%s''.', field);
end
if(~isempty(err)), throw(err); end

function printfields()
%Print out fields with defaults

fprintf('  CBC GENERAL SETTINGS:\n');
fprintf('                  dualTol: [ Dual Infeasibility Tolerance: {1e-7} ]\n');
fprintf('                primalTol: [ Primal Infeasibility Tolerance: {1e-7} ]\n');
fprintf('             allowableGap: [ Allowable gap between best known solution and best possible: {0} ]\n');
fprintf('         allowableFracGap: [ Allowable fraction gap between best known solution and best possible: {0} ]\n');

fprintf('\n  CBC MIP SETTINGS:\n');
fprintf('                 startalg: [ Root Node algorithm: Primal Simplex (''Primal''), Dual Simplex (''Dual''), Barrier (''Barrier''), Default {} ]\n');
fprintf('                 strategy: [ Switches on groups of features: Easy problems (0), Default {1}, Aggressive (2) ]\n');
fprintf('               preprocess: [ Whether to use integer preprocessing: ''Off'', ''On'', ''Equal'', {''SOS''}, ''TrySOS'', ''EqualAll'', ''Strategy'', ''Aggregate'', ''ForceSOS'' ]\n');
fprintf('               heuristics: [ Switches most heuristics: ''Off'', {''On''} ]\n');
fprintf('             nodeStrategy: [ What strategy to use to select nodes: ''Hybrid'', {''Fewest''}, ''Depth'', ''UpFewest'', ''DownFewest'', ''UpDepth'', ''DownDepth'' ]\n');
fprintf('                increment: [ A valid solution must be at least this much better than the last integer solution: {1e-5} ]\n');
fprintf('             maxSolutions: [ Maximum number of solutions to get: Default {-1}, User (> 0) ]\n');
fprintf('          strongBranching: [ Maximum number of candidates to be evaluated for strong branching: {5} ]\n');
fprintf('         trustPseudoCosts: [ Number of branches before we trust pseudocosts: {5} ]\n');
fprintf('          combineSolution: [ Whether to use combine solution heuristic: ''Off'', {''On''} ]\n');
fprintf('                     dins: [ Whether to try Distance Induced Neighbourhood Search: {''Off''}, ''On'' ]\n');
fprintf('               divingSome: [ Whether to try Diving heuristics: {''Off''}, ''On'', ''Both'', ''Before'' ]\n');
fprintf('        divingCoefficient: [ Whether to try DiveCoefficient: ''Off'', {''On''}, ''Both'', ''Before'' ]\n');
fprintf('         divingFractional: [ Whether to try DiveFractional: {''Off''}, ''On'', ''Both'', ''Before'' ]\n');
fprintf('             divingGuided: [ Whether to try DiveGuided: {''Off''}, ''On'', ''Both'', ''Before'' ]\n');
fprintf('         divingLineSearch: [ Whether to try DiveLineSearch: {''Off''}, ''On'', ''Both'', ''Before'' ]\n');
fprintf('         divingPseudoCost: [ Whether to try DivePseudoCost: {''Off''}, ''On'', ''Both'', ''Before'' ]\n');
fprintf('       divingVectorLength: [ Whether to try DiveVectorLength: {''Off''}, ''On'', ''Both'', ''Before'' ]\n');
fprintf('          feasibilityPump: [ Whether to try Feasibility Pump: ''Off'', {''On''}, ''Both'', ''Before'' ]\n');
fprintf('    feasibilityPumpPasses: [ How many passes in Feasibility Pump: {30} ]\n');
fprintf('          greedyHeuristic: [ Whether to use a greedy heuristic: {''Off''}, ''On'', ''Both'', ''Before'' ]\n');
fprintf('          localTreeSearch: [ Whether to use local treesearch: {''Off''}, ''On'' ]\n');
fprintf('              pivotAndFix: [ Whether to try Pivot and Fix heuristic: {''Off''}, ''On'', ''Both'', ''Before'' ]\n');
fprintf('       randomizedRounding: [ Whether to try Randomized Rounding heuristic: {''Off''}, ''On'', ''Both'', ''Before'' ]\n');
fprintf('                     rens: [ Whether to try Relaxation Enforced Neighbourhood Search: {''Off''}, ''On'', ''Both'', ''Before'', ''200'', ''1000'', ''10000'', ''dj'', ''djbefore'' ]\n');
fprintf('                     rins: [ Whether to try Relaxed Induced Neighbourhood Search: ''Off'', {''On''}, ''Both'', ''Before'' ]\n');
fprintf('        roundingHeuristic: [ Whether to try Rounding heuristic: ''Off'', {''On''}, ''Both'', ''Before'' ]\n');
fprintf('             vubHeuristic: [ Type of vub heuristic: {-1} ]\n');
fprintf('             costStrategy: [ How to use costs as priorities: {''Off''}, ''Priorities'', ''ColumnOrder'', ''01first'', ''01last'', ''Length'', ''Singletons'', ''Nonzero'' ]\n');

fprintf('\n  CBC CUT SETTINGS:\n');
fprintf('                  useCuts: [ Switches all cuts on or off: ''Off'', {''On''}, ''Root'', ''IfMove'', ''ForceOn'' ]\n');
fprintf('                 cutDepth: [ Depth in tree at which to do cuts: Default {-1}, User (>= 0) ]\n');
fprintf('         maxCutPassesRoot: [ Maximum number of cut passes at the root node: Default {-1}, User (>= 0) ]\n');
fprintf('             maxCutPasses: [ Maximum number of cut passes at other nodes: {1} ]\n');
fprintf('            slowCutPasses: [ Maximum number of tries for slower cuts: {10} ]\n');

fprintf('               cliqueCuts: [ Whether to use Clique cuts: ''Off'', ''On'', ''Root'', {''IfMove''}, ''ForceOn'', ''OnGlobal'' ]\n');
fprintf('            flowCoverCuts: [ Whether to use Flow Cover cuts: ''Off'', ''On'', ''Root'', {''IfMove''}, ''ForceOn'', ''OnGlobal'' ]\n');
fprintf('                  GMICuts: [ Whether to use alternative Gomory cuts: {''Off''}, ''On'', ''Root'', ''IfMove'', ''ForceOn'', ''EndOnly'', ''Long'', ''LongRoot'', ''LongIfMove'', ''ForceLongOn'', ''LongEndOnly'' ]\n');
fprintf('             lagomoryCuts: [ Whether to use Lagrangean Gomory cuts: {''Off''}, ''EndOnlyRoot'', ''EndCleanRoots'', ''EndBothRoot'', ''EndOnly'', ''EndClean'', ''EndBoth'', ''OnlyAsWell'', ''CleanAsWell'', ''BothAsWell'', ''OnlyInstead'', ''CleanInstead'', ''BothInstead'' ]\n');
fprintf('               gomoryCuts: [ Whether to use Gomory cuts: ''Off'', ''On'', ''Root'', {''IfMove''}, ''ForceOn'', ''OnGlobal'', ''ForceAndGlobal'', ''ForceLongOn'', ''Long'' ]\n');
fprintf('             knapsackCuts: [ Whether to use Knapsack cuts: ''Off'', ''On'', ''Root'', {''IfMove''}, ''ForceOn'', ''OnGlobal'', ''ForceAndGlobal'' ]\n');
fprintf('       liftAndProjectCuts: [ Whether to use Lift and Project cuts: {''Off''}, ''On'', ''Root'', ''IfMove'', ''ForceOn'' ]\n');
fprintf(' mixedIntegerRoundingCuts: [ Whether to use Mixed Integer Rounding cuts: ''Off'', ''On'', ''Root'', {''IfMove''}, ''ForceOn'', ''OnGlobal'' ]\n');
fprintf('              probingCuts: [ Whether to use Probing cuts: ''Off'', {''On''}, ''Root'', ''IfMove'', ''ForceOn'', ''OnGlobal'', ''ForceOnGlobal'', ''ForceOnBut'', ''ForceOnStrong'', ''ForceOnButStrong'', ''StrongRoot'' ]\n');
fprintf('       reduceAndSplitCuts: [ Whether to use Reduce-and-Split cuts: {''Off''}, ''On'', ''Root'', ''IfMove'', ''ForceOn'' ]\n');
fprintf('     residualCapacityCuts: [ Whether to use Residual Capacity cuts: {''Off''}, ''On'', ''Root'', ''IfMove'', ''ForceOn'' ]\n');
fprintf('               twoMirCuts: [ Whether to use Two Phase Mixed Integer Rounding cuts: ''Off'', ''On'', {''Root''}, ''IfMove'', ''ForceOn'', ''OnGlobal'', ''ForceAndGlobal'', ''ForceLongOn'' ]\n');
fprintf('             latwoMirCuts: [ Whether to use Lagrangean TwoMir cuts: {''Off''}, ''EndOnlyRoot'', ''EndCleanRoot'', ''EndBothRoot'', ''EndOnly'', ''EndClean'', ''EndBoth'', ''OnlyAsWell'', ''CleanAsWell'', ''BothAsWell'', ''OnlyInstead'', ''CleanInstead'', ''BothInstead'' ]\n');
fprintf('             zeroHalfCuts: [ Whether to use Zero Half cuts: {''Off''}, ''On'', ''Root'', ''IfMove'', ''ForceOn'', ''OnGlobal'' ]\n');

fprintf('\n  CLP RELAXED SOLVER SETTINGS:\n');
fprintf('               idiotCrash: [ Whether to try idiot crash: Let Code Decide {-1}, Off (0), n passes (>0) ]\n');
fprintf('              sprintCrash: [ Whether to try sprint crash: Let Code Decide {-1}, Off (0), n passes (>0) ]\n');
fprintf('                    crash: [ Whether to create basis for problem: {''Off''}, ''On'', ''solow_halim'', ''lots'' ]\n');
fprintf('            factorization: [ Factorization Selection: {''Normal''}, ''Dense'', ''Simple'', ''OSL'' ]\n');
fprintf('                maxFactor: [ Maximum iterations between refactorizations: Clp heuristic {200}, User Set (>0) ]\n');
fprintf('                crossover: [ Whether to get a basic solution after barrier: {''On''}, ''Off'', ''Maybe'', ''Presolve'' ]\n');
fprintf('                dualPivot: [ Dual pivot choice algorithm: {''automatic''}, ''dantzig'', ''partial'', ''steepest'' ]\n');
fprintf('              primalPivot: [ Primal pivot choice algorithm: {''automatic''}, ''exact'', ''dantzig'', ''partial'', ''steepest'', ''change'', ''sprint'' ]\n');
fprintf('                dualBound: [ Initially algorithm acts as if no gap between bounds exceeds this value: {1e10} ]\n');
fprintf('             primalWeight: [ Initially algorithm acts as if it costs this much to be infeasible: {1e10} ]\n');
fprintf('             perturbation: [ Whether to perturb problem: {''On''}, ''Off'' ]\n');
fprintf('                  scaling: [ Whether to scale the problem: {''automatic''}, ''equilibrium'', ''geometric'', ''dynamic'', ''rowsonly'' ]\n');
fprintf('                 presolve: [ Whether to presolve the problem: {''On''}, ''Off'', ''More'' ]\n');
fprintf('             preTolerance: [ Tolerance to use in presolve: {1e-8} ]\n');
fprintf('             passPresolve: [ How many passes in presolve: {5} ]\n');

fprintf('\n');
