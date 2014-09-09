function options = bonminset(varargin)
%BONMINSET  Create or alter the options for Optimization with BONMIN
%
% options = bonminset('param1',value1,'param2',value2,...) creates an IPOPT
% options structure with the parameters 'param' set to their corresponding
% values in 'value'. Parameters not specified will be set to the IPOPT
% default.
%
% options = bonminset(oldopts,'param1',value1,...) creates a copy of the old
% options 'oldopts' and then fills in (or writes over) the parameters
% specified by 'param' and 'value'.
%
% options = bonminset() creates an options structure with all fields set to
% bonminset defaults.
%
% bonminset() prints a list of all possible fields and their function.
%
% See supplied BONMIN and IPOPT Documentation for further details of these options.

%   Copyright (C) 2011 Jonathan Currie (I2C2)

% Print out possible values of properties.
if((nargin == 0) && (nargout == 0))
    printfields();
    return
end

%Names and Defaults
Names = {'algorithm','var_lin','cons_lin','allowable_fraction_gap','allowable_gap','cutoff','cutoff_decr','enable_dynamic_nlp','iteration_limit',...
         'nlp_failure_behavior','node_comparison','num_cut_passes','num_cut_passes_at_root',...
         'number_before_trust','number_strong_branch','solution_limit','sos_constraints','tree_search_strategy',...
         'variable_selection','feasibility_pump_objective_norm','heuristic_RINS','heuristic_dive_MIP_fractional',...
         'heuristic_dive_MIP_vectorLength','heuristic_dive_fractional','heuristic_dive_vectorLength',...
         'heuristic_feasibility_pump','pump_for_minlp','candidate_sort_criterion','maxmin_crit_have_sol',...
         'maxmin_crit_no_sol','min_number_strong_branch','number_before_trust_list','number_look_ahead',...
         'number_strong_branch_root','setup_pseudo_frac','coeff_var_threshold','dynamic_def_cutoff_decr','first_perc_for_cutoff_decr',...
         'max_consecutive_infeasible','num_resolve_at_infeasibles',...
         'num_resolve_at_node','num_resolve_at_root','second_perc_for_cutoff_decr','milp_solver','milp_strategy','bb_log_interval',...
         'bb_log_level','lp_log_level','milp_log_level','nlp_log_level','nlp_log_at_root','oa_log_frequency','oa_log_level',...
         'max_consecutive_failures','max_random_point_radius','num_retry_unsolved_random_point','random_point_perturbation_interval',...
         'random_point_type','warm_start'}';
Defaults = {'B-BB',[],[],0,0,1e100,1e-5,'no',2^31-1,'stop','best-bound',1,20,8,20,2^31-1,'enable','probed-dive','osi-strong',1,'no',...
            'yes','no','no','no','no','no','best-ps-cost',0.1,0.7,0,0,0,2^31-1,0.5,0.1,'no',-0.02,0,0,0,0,-0.05,'Cbc_D','find_good_sol',...
            100,1,0,0,0,0,100,1,10,100000,0,1,'Jon','none'}';  
%Get IPOPT Defaults if required
if(nargin == 2 && strcmpi(varargin{2},'noipopt'))
    varargin = varargin(1);
else
    iset = ipoptset('bonmin');       
    Names = [Names;fieldnames(iset)];
    Defaults = [Defaults;struct2cell(iset)];
end

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
    case {'first_perc_for_cutoff_decr','second_perc_for_cutoff_decr'}
        err = opticheckval.checkScalarDbl(value,field);
    %Scalar non negative double
    case {'min_number_strong_branch','number_look_ahead','number_strong_branch_root','coeff_var_threshold'}
        err = opticheckval.checkScalarNonNeg(value,field);      
    %Scalar non zero double
    case {'oa_log_frequency','max_random_point_radius','random_point_perturbation_interval'}
        err = opticheckval.checkScalarGrtZ(value,field);    
    %Scalar non negative integer
    case {'iteration_limit','num_cut_passes','num_cut_passes_at_root','number_before_trust','number_strong_branch',...
          'solution_limit','max_consecutive_infeasible','num_resolve_at_infeasibles','num_resolve_at_node',...
          'num_resolve_at_root','bb_log_interval','max_consecutive_failures','num_retry_unsolved_random_point'}
        err = opticheckval.checkScalarIntNonNeg(value,field);      
    %Scalar double with bounds
    case 'cutoff_decr'
        err = opticheckval.checkScalarBoundLELE(value,field,-1e10,1e10);
    case {'allowable_fraction_gap','allowable_gap'}
        err = opticheckval.checkScalarBoundLELE(value,field,-1e20,1e20);
    case 'cutoff'
        err = opticheckval.checkScalarBoundLELE(value,field,-1e100,1e100);
    case {'maxmin_crit_have_sol','maxmin_crit_no_sol','setup_pseudo_frac'}
        err = opticheckval.checkScalarBoundLELE(value,field,0,1);
    case 'number_before_trust_list'
        err = opticheckval.checkScalarBoundLEL(value,field,-1,inf);
    %Scalar integer with bounds
    case 'feasibility_pump_objective_norm'
        err = opticheckval.checkScalarIntBoundLELE(value,field,1,2);
    case 'bb_log_level'
        err = opticheckval.checkScalarIntBoundLELE(value,field,0,5);
    case {'lp_log_level','milp_log_level'}
        err = opticheckval.checkScalarIntBoundLELE(value,field,0,4);
    case {'nlp_log_level','oa_log_level'}
        err = opticheckval.checkScalarIntBoundLELE(value,field,0,2);
    case 'nlp_log_at_root'
        err = opticheckval.checkScalarIntBoundLELE(value,field,0,12);
    %Vector integer 0/1
    case {'var_lin','cons_lin'}
        err = opticheckval.checkVector01(value,field);
    %Yes / No
    case {'heuristic_rins','heuristic_dive_mip_fractional','heuristic_dive_mip_vectorlength','heuristic_dive_fractional',...
          'heuristic_dive_vectorlength','heuristic_feasibility_pump','pump_for_minlp','dynamic_def_cutoff_decr'}
        err = opticheckval.checkYesNo(value,field);
    %Yes / No with Yes Crash
    case 'enable_dynamic_nlp'
        err = opticheckval.checkYesCrashNo(value,field);
    %Enable / Disable
    case 'sos_constraints'
        err = opticheckval.checkEnDis(value,field);        
    %Misc String methods
    case 'algorithm'
        err = opticheckval.checkValidString(value,field,{'B-BB','B-OA','B-QG','B-Hyb','B-Ecp','B-iFP'});
    case 'nlp_failure_behavior'
        err = opticheckval.checkValidString(value,field,{'stop','fathom'});
    case 'node_comparison'
        err = opticheckval.checkValidString(value,field,{'best-bound','depth-first','breadth-first','dynamic','best-guess'});
    case 'tree_search_strategy'
        err = opticheckval.checkValidString(value,field,{'top-node','dive','probed-dive','dfs-dive','dfs-dive-dynamic'});
    case 'variable_selection'
        err = opticheckval.checkValidString(value,field,{'most-fractional','strong-branching','reliability-branching','qp-strong-branching',...
                                                         'lp-strong-branching','nlp-strong-branching','osi-simple','osi-strong','random'});
    case 'candidate_sort_criterion'
        err = opticheckval.checkValidString(value,field,{'best-ps-cost','worst-ps-cost','most-fractional','least-fractional'});    
    case 'milp_solver'
        err = opticheckval.checkValidString(value,field,{'cbc_d','cplex'});   
    case 'milp_strategy'
        err = opticheckval.checkValidString(value,field,{'find_good_sol','solve_to_optimality'}); 
    case 'random_point_type'
        err = opticheckval.checkValidString(value,field,{'jon','andreas','claudia'});
    case 'warm_start'
        err = opticheckval.checkValidString(value,field,{'none','fake_basis','optimum','interior_point'});        
        
    %IPOPT or Unknown Parameter
    otherwise  
        err = [];
        try
            %try for ipoptset arg
            iset = ipoptset(field,value);            
        catch ME
            throw(ME);
        end
        if(isempty(iset))
            err = MException('OPTI:SetFieldError','Unrecognized parameter name ''%s''.', field);
        end
end
if(~isempty(err)), throw(err); end


function printfields()
%Print out fields with defaults
fprintf(' ----- BONMIN (MIP SOLVER) SETTINGS ----- \n');
fprintf('  BONMIN GENERAL SETTINGS:\n');
fprintf('                          algorithm: [ MIP Algorithm {''B-BB''}, ''B-OA'', ''B-QG'', ''B-Hyb'', ''B-Ecp'', ''B-iFP'' ] \n');
fprintf('                            var_lin: [ Decision Variable Linearity (0 - Nonlinear, 1 - Linear) {zeros(n,1)} ] \n');
fprintf('                           cons_lin: [ Constraint Linearity (0 - Nonlinear, 1 - Linear) {zeros(m,1)} ] \n');

fprintf('\n  BONMIN BRANCH & BOUND SETTINGS:\n');
fprintf('             allowable_fraction_gap: [ Value of Relative Gap when algorithm stops {0} ] \n');
fprintf('                      allowable_gap: [ Value of Absolute Gap when algorithm stops {0} ] \n');
fprintf('                             cutoff: [ The algorithm will only look for values better than cutoff {1e100} ] \n');
fprintf('                        cutoff_decr: [ Amount by which cutoff is decremented below a new best upper bound {1e-5} ] \n');
fprintf('                 enable_dynamic_nlp: [ Enable dynamic linear and quadratic rows addition in NLP ''yes'', {''no''} ] \n');
fprintf('                    iteration_limit: [ Cumulated maximum number of iterations in the algorithm to process continuous relaxations {2^31-1} ] \n');
fprintf('               nlp_failure_behavior: [ Action to take when an NLP is unsolved by IPOPT {''stop''}, ''fathom'' ] \n');
fprintf('                    node_comparison: [ Node selection strategy {''best-bound''}, ''depth-first'', ''breadth-first'', ''dynamic'', ''best-guess'' ] \n');
fprintf('                     num_cut_passes: [ Maximum number of cut passes at regular nodes {1} ] \n');
fprintf('             num_cut_passes_at_root: [ Maximum number of cut passes at root node {20} ] \n');
fprintf('                number_before_trust: [ Number of branches on a variable before its pseudo costs are to be believed in dynamic strong branching {8} ] \n');
fprintf('               number_strong_branch: [ Maximum number of variables considered for strong branching {20} ] \n');
fprintf('                     solution_limit: [ Abort after this many integer feasible solutions have been found (or 0 to deactivate) {2^31-1} ] \n');
fprintf('                    sos_constraints: [ Activate SOS Type 1 constraints {''enable''}, ''disable'' ] \n');
fprintf('               tree_search_strategy: [ Tree traversing strategy {''probed-dive''}, ''top-node'', ''dive'', ''dfs-dive'', ''dfs-dive-dynamic'' ] \n');
fprintf('                 variable_selection: [ Variable selection strategy {''osi-strong''}, ''most-fractional'', ''reliability-branching'', ''qp-strong-branching'', ''lp-strong-branching'', ''nlp-strong-branching'', ''osi-simple'', ''strong-branching'', ''random'' ] \n');

fprintf('\n  BONMIN MINLP HEURISTIC SETTINGS:\n');
fprintf('    feasibility_pump_objective_norm: [ Norm of feasibility pump objective norm {1}, 2 ] \n');
fprintf('                     heuristic_RINS: [ Use RINS heuristic {''no''}, ''yes'' ] \n');
fprintf('      heuristic_dive_MIP_fractional: [ Use Dive MIP fractional heuristic {''no''}, ''yes'' ] \n');
fprintf('    heuristic_dive_MIP_vectorLength: [ Use Dive MIP vectorLength heuristic {''no''}, ''yes'' ] \n');
fprintf('          heuristic_dive_fractional: [ Use Dive fractional heuristic ''no'', {''yes''} ] \n');
fprintf('        heuristic_dive_vectorLength: [ Use Dive vectorLength heuristic ''no'', {''yes''} ] \n');
fprintf('         heuristic_feasibility_pump: [ Use feasibility pump heuristic {''no''}, ''yes'' ] \n');
fprintf('                     pump_for_minlp: [ Use FP for MINLP {''no''}, ''yes'' ] \n');

fprintf('\n  BONMIN NON-CONVEX PROBLEM SETTINGS:\n');
fprintf('                coeff_var_threshold: [ Coefficient of variation threshold (for dynamic definition of cutoff_decr) {0.1} ] \n');
fprintf('            dynamic_def_cutoff_decr: [ Define cutoff_decr dynamically? {''no''}, ''yes'' ] \n');
fprintf('         first_perc_for_cutoff_decr: [ Percentage used when the coeff of variance is smaller than the threshold, to compute cutoff_decr dynamically {-0.02} ] \n');
fprintf('         max_consecutive_infeasible: [ Maximum number of consecutive infeasible subproblems before aborting a branch {0} ] \n');
fprintf('         num_resolve_at_infeasibles: [ Number of tries to resolve an infeasible node (other than the root) of the tree with different starting points {0} ] \n');
fprintf('                num_resolve_at_node: [ Number of tries to resolve a node (other than the root) of the tree with different starting points (like multi-start) {0} ] \n');
fprintf('                num_resolve_at_root: [ Number of tries to resolve the root node with different starting points {0} ] \n');
fprintf('        second_perc_for_cutoff_decr: [ Percentage used when the coeff of variance is greater than the threshold, to compute cutoff_decr dynamically {-0.05} ] \n');

fprintf('\n  BONMIN NLP SOLUTION ROBUSTNESS SETTINGS:\n');
fprintf('           max_consecutive_failures: [ Number of consecutive unsolved problems before aborting a branch of the tree {10} ] \n');
fprintf('            max_random_point_radius: [ Max value for coordinate of a random point {100000} ] \n');
fprintf('    num_retry_unsolved_random_point: [ Number of times that the algorithm will try to resolve an unsolved NLP with a random starting point (with new point) {0} ] \n');
fprintf(' random_point_perturbation_interval: [ Amount by which starting point is perturbed when choosing a random point {1} ] \n');
fprintf('                  random_point_type: [ Random starting point method {''Jon''}, ''Andreas'', ''Claudia'' ] \n');
fprintf('                         warm_start: [ Warm start method {''none''}, ''fake_basis'', ''optimum'', ''interior'' ] \n');

fprintf('\n  BONMIN STRONG BRANCHING SETTINGS:\n');
fprintf('           candidate_sort_criterion: [ Criterion used to choose candidates in strong-branching {''best-ps-cost''}, ''worst-ps-cost'', ''most-fractional'', ''least-fractional'' ] \n')
fprintf('               maxmin_crit_have_sol: [ Weight towards minimum of lower and upper branching estimates when a solution has been found {0.1} ] \n')
fprintf('                 maxmin_crit_no_sol: [ Weight towards minimum of lower and upper branching estimates when no solution has been found {0.7} ] \n')
fprintf('           min_number_strong_branch: [ Minimum number of variables for strong branching (overriding trust) {0} ] \n')
fprintf('           number_before_trust_list: [ Number of branches on a variable vbefore its pseudo costs are to be believed during setup of strong branching list {0} ] \n')
fprintf('                  number_look_ahead: [ Limit of look-ahead strong-branching trials {0} ] \n')
fprintf('          number_strong_branch_root: [ Maximum number of variables considered for strong branching at root node {2^31-1} ] \n')
fprintf('                  setup_pseudo_frac: [ Proportion of strong branching list that has to be taken from most integer infeasible list {0.5} ] \n')

fprintf('\n  BONMIN MILP SOLVER SETTINGS:\n');
fprintf('                        milp_solver: [ Subsolver to solve MILP problems in OA decompositions {''Cbc_D''}, ''Cplex'' (must be installed on your PC) ] \n');
fprintf('                      milp_strategy: [ MILP solving strategy {''find_good_sol''}, ''solve_to_optimality'' ] \n');

fprintf('\n  BONMIN OUPUT AND LOG LEVEL SETTINGS:\n');
fprintf('                    bb_log_interval: [ Interval at which node level output is printed (number of nodes) {100} ] \n');
fprintf('                       bb_log_level: [ Level of branch and bound log detail (0-5) {1} ] \n');
fprintf('                       lp_log_level: [ Level of LP solver log detail (0-4) {0} ] \n');
fprintf('                     milp_log_level: [ Level of MILP solver log detail (0-4) {0} ] \n');
fprintf('                      nlp_log_level: [ Level of NLP solver log detail (0-2) {0} ] \n');
fprintf('                    nlp_log_at_root: [ Level of NLP solver log detail at root node (0-12) {0} ] \n');
fprintf('                   oa_log_frequency: [ Frequency (in seconds) of OA log messages  {100} ] \n');
fprintf('                       oa_log_level: [ Level of OA decomposition log detail (0-2) {1} ] \n');

fprintf('\n ----- IPOPT (RELAXED SOLVER) SETTINGS ----- \n');
ipoptset('bonmin');

