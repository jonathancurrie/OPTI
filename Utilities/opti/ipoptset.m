function options = ipoptset(varargin)
%IPOPTSET  Create or alter the options for Optimization with IPOPT
%
% options = ipoptset('param1',value1,'param2',value2,...) creates an IPOPT
% options structure with the parameters 'param' set to their corresponding
% values in 'value'. Parameters not specified will be set to the IPOPT
% default.
%
% options = ipoptset(oldopts,'param1',value1,...) creates a copy of the old
% options 'oldopts' and then fills in (or writes over) the parameters
% specified by 'param' and 'value'.
%
% options = ipoptset() creates an options structure with all fields set to
% IPOPTSET defaults.
%
% ipoptset() prints a list of all possible fields and their function.
%
% See supplied IPOPT Documentation for further details of these options.

%   Copyright (C) 2011-2013 Jonathan Currie (Control Engineering)

%Default mode (bonmin mode uses slightly different options)
mode = 'ipopt';

% Print out possible values of properties.
if((nargin == 0) && (nargout == 0))
    printfields('ipopt');
    return
elseif(nargin == 1) && (strcmpi(varargin{1},'bonmin'))    
    if(nargout == 0)
        printfields('bonmin');
        return;        
    else
        mode = 'bonmin';
        varargin = {};
    end
end
%Names and Defaults
Names = {'dual_inf_tol','constr_viol_tol','compl_inf_tol','acceptable_tol','acceptable_iter','acceptable_constr_viol_tol',...
        'acceptable_dual_inf_tol','acceptable_compl_inf_tol','acceptable_obj_change_tol','diverging_iterates_tol',...
        'obj_scaling_factor','nlp_scaling_method','nlp_scaling_max_gradient','nlp_scaling_min_value',...
        'bound_relax_factor','honor_original_bounds','check_derivatives_for_naninf','nlp_lower_bound_inf',...
        'nlp_upper_bound_inf','fixed_variable_treatment','jac_c_constant','jac_d_constant','hessian_constant',...
        'bound_frac','bound_push','slack_bound_frac','slack_bound_push','bound_mult_init_val','constr_mult_init_max',...
        'bound_mult_init_method','mehrotra_algorithm','mu_strategy','mu_oracle','quality_function_max_section_steps',...
        'fixed_mu_oracle','adaptive_mu_globalization','mu_init','mu_max_fact','mu_max','mu_min','mu_target',...
        'barrier_tol_factor','mu_linear_decrease_factor','mu_superlinear_decrease_power'...
        'alpha_for_y','alpha_for_y_tol','recalc_y','recalc_y_feas_tol','max_soc','watchdog_shortened_iter_trigger','watchdog_trial_iter_max',...
        'accept_every_trial_step','warm_start_init_point','warm_start_bound_push','warm_start_bound_frac','warm_start_slack_bound_push',...
        'warm_start_slack_bound_frac','warm_start_mult_bound_push','warm_start_mult_init_max','expect_infeasible_problem',...
        'expect_infeasible_problem_ctol','expect_infeasible_problem_ytol','start_with_resto','soft_resto_pderror_reduction_factor',...
        'required_infeasibility_reduction','bound_mult_reset_threshold','constr_mult_reset_threshold','evaluate_orig_obj_at_resto_trial',...
        'linear_solver','linear_system_scaling','linear_scaling_on_demand','max_refinement_steps','min_refinement_steps',...
        'max_hessian_perturbation','min_hessian_perturbation','first_hessian_perturbation','perturb_inc_fact_first','perturb_inc_fact',...
        'perturb_dec_fact','jacobian_regularization_value','hessian_approximation','limited_memory_max_history',...
        'limited_memory_max_skipping','limited_memory_initialization','limited_memory_init_val','limited_memory_init_val_max',...
        'limited_memory_init_val_min','limited_memory_special_for_resto','derivative_test','derivative_test_perturbation',...
        'derivative_test_tol','derivative_test_print_all','derivative_test_first_index','point_perturbation_radius',...
        'ma57_pivtol','ma57_pivtolmax','ma57_pre_alloc','ma57_pivot_order','ma57_automatic_scaling','ma57_block_size',...
        'ma57_node_amalgamation','ma57_small_pivot_flag','mumps_pivtol','mumps_pivtolmax','mumps_mem_percent','mumps_permuting_scaling',...
        'mumps_pivot_order','mumps_scaling','pardiso_redo_symbolic_fact_only_if_inertia_wrong','pardiso_order'};

Defaults = {1,0.0001,0.0001,1e-6,15,0.01,1e10,0.01,1e20,1e20,1,'gradient-based',100,1e-8,1e-8,'yes','no',-1e19,1e19,...
        'make_parameter','no','no','no',0.01,0.01,0.01,0.01,1,1000,'constant','no','monotone','quality-function',...
        8,'average_compl','obj-constr-filter',0.1,1000,100000,1e-11,0,10,0.2,1.5,'primal',10,'no',1e-6,...
        4,10,3,'no','no',0.001,0.001,0.001,0.001,0.001,1e-6,'no',0.001,1e8,'no',0.9999,0.9,1000,0,'yes',...
        'ma57','none','yes',10,1,1e20,1e-20,0.0001,100,8,1/3,1e-8,'exact',6,2,'scalar1',1,1e8,1e-8,'no',...
        'none',1e-8,0.0001,'no',-2,10,1e-8,0.0001,1.05,5,'yes',16,16,0,1e-6,0.1,1000,7,7,77,'no','metis'};

%Modify defaults if setting up for BONMIN
if(strcmpi(mode,'bonmin'))
    ind = strcmp(Defaults,'quality-function');
    Defaults{ind} = 'probing';
    ind = strcmp(Names,'required_infeasibility_reduction');
    Defaults{ind} = 0.1;
    ind = strcmp(Names,'expect_infeasible_problem');
    Defaults{ind} = 'yes';
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
    case {'obj_scaling_factor','nlp_lower_bound_inf','nlp_upper_bound_inf','warm_start_mult_init_max'}
        err = opticheckval.checkScalarDbl(value,field);
    %Scalar non negative double
    case {'acceptable_obj_change_tol','nlp_scaling_min_value','bound_relax_factor','constr_mult_init_max',...
          'quality_function_max_section_steps','mu_target','alpha_for_y_tol','expect_infeasible_problem_ctol',...
          'soft_resto_pderror_reduction_factor','bound_mult_reset_threshold','constr_mult_reset_threshold',...
          'min_hessian_perturbation','jacobian_regularization_value','point_perturbation_radius','warm_start_mult_bound_push'}
        err = opticheckval.checkScalarNonNeg(value,field);      
    %Scalar non zero double
    case {'dual_inf_tol','constr_viol_tol','compl_inf_tol','acceptable_tol','acceptable_constr_viol_tol',...
          'acceptable_dual_inf_tol','acceptable_compl_inf_tol','diverging_iterates_tol','nlp_scaling_max_gradient',...
          'bound_push','slack_bound_push','bound_mult_init_val','mu_init','mu_max_fact','mu_max','mu_min',...
          'barrier_tol_factor','recalc_y_feas_tol','warm_start_bound_push','warm_start_slack_bound_push',...
          'expect_infeasible_problem_ytol','max_hessian_perturbation','first_hessian_perturbation',...
          'limited_memory_init_val','limited_memory_init_val_max','limited_memory_init_val_min',...
          'derivative_test_perturbation','derivative_test_tol'}
        err = opticheckval.checkScalarGrtZ(value,field);    
    %Scalar non negative integer
    case {'acceptable_iter','max_soc','watchdog_shortened_iter_trigger','limited_memory_max_history',...
          'mumps_mem_percent','max_refinement_steps','min_refinement_steps'}
        err = opticheckval.checkScalarIntNonNeg(value,field);  
    %Scalar double with bounds
    case {'bound_frac','slack_bound_frac','warm_start_bound_frac','warm_start_slack_bound_frac'}
        err = opticheckval.checkScalarBoundLLE(value,field,0,0.5);
    case 'required_infeasibility_reduction'
        err = opticheckval.checkScalarBoundLEL(value,field,0,1);
    case 'derivative_test_first_index'
        err = opticheckval.checkScalarBoundLEL(value,field,-2,inf);
    case 'ma57_pre_alloc'
        err = opticheckval.checkScalarBoundLEL(value,field,1,inf);
    case {'mu_linear_decrease_factor','perturb_dec_fact','ma57_pivtol','ma57_pivtolmax','mumps_pivtol',...
          'mumps_pivtolmax'}
        err = opticheckval.checkScalarBoundLL(value,field,0,1);
    case 'mu_superlinear_decrease_power'
        err = opticheckval.checkScalarBoundLL(value,field,1,2);
    case {'perturb_inc_fact_first','perturb_inc_fact'}
        err = opticheckval.checkScalarBoundLL(value,field,1,inf);
    %Scalar integer with bounds
    case {'watchdog_trial_iter_max','limited_memory_max_skipping','ma57_block_size','ma57_node_amalgamation'}
        err = opticheckval.checkScalarIntBoundLEL(value,field,1,inf);
    case 'ma57_pivot_order'
        err = opticheckval.checkScalarIntBoundLELE(value,field,0,5);
    case 'ma57_small_pivot_flag'
        err = opticheckval.checkScalarIntBoundLELE(value,field,0,1);
    case {'mumps_permuting_scaling','mumps_pivot_order'}
        err = opticheckval.checkScalarIntBoundLELE(value,field,0,7);
    case 'mumps_scaling'
        err = opticheckval.checkScalarIntBoundLELE(value,field,-2,77);
    %Yes / No
    case {'honor_original_bounds','check_derivatives_for_naninf','jac_c_constant','jac_d_constant',...
          'hessian_constant','mehrotra_algorithm','recalc_y','accept_every_trial_step','warm_start_init_point',...
          'expect_infeasible_problem','start_with_resto','evaluate_orig_obj_at_resto_trial','linear_scaling_on_demand',...
          'limited_memory_special_for_resto','derivative_test_print_all','ma57_automatic_scaling','pardiso_redo_symbolic_fact_only_if_inertia_wrong'}
        err = opticheckval.checkYesNo(value,field);
    %Misc String methods
    case 'nlp_scaling_method'
        err = opticheckval.checkValidString(value,field,{'none','gradient-based'});
    case 'fixed_variable_treatment'
        err = opticheckval.checkValidString(value,field,{'make_parameter','make_constraint','relax_bounds'});
    case 'bound_mult_init_method'
        err = opticheckval.checkValidString(value,field,{'constant','mu-based'});
    case 'mu_strategy'
        err = opticheckval.checkValidString(value,field,{'monotone','adaptive'});
    case 'mu_oracle'
        err = opticheckval.checkValidString(value,field,{'probing','loqo','quality-function'});
    case 'fixed_mu_oracle'
        err = opticheckval.checkValidString(value,field,{'probing','loqo','quality-function','average_compl'});
    case 'adaptive_mu_globalization'
        err = opticheckval.checkValidString(value,field,{'kkt-error','obj-constr-filter','never-monotone-mode'});
    case 'alpha_for_y'
        err = opticheckval.checkValidString(value,field,{'primal','bound-mult','min','max','full','min-dual-infeas','safer-min-dual-infeas','primal-and-full','dual-and-full','acceptor'});
    case 'linear_solver'
        err = opticheckval.checkValidString(value,field,{'ma57','pardiso','mumps'});
    case 'linear_system_scaling'
        err = opticheckval.checkValidString(value,field,{'none','slack-based'});
    case 'hessian_approximation'
        err = opticheckval.checkValidString(value,field,{'exact','limited-memory'});
    case 'limited_memory_initialization'
        err = opticheckval.checkValidString(value,field,{'scalar1','scalar2','scalar3','scalar4','constant'});
    case 'derivative_test'
        err = opticheckval.checkValidString(value,field,{'none','first-order','second-order','only-second-order'});
    case 'pardiso_order'
        err = opticheckval.checkValidString(value,field,{'metis','amd','pmetis'});
    %Unknown Parameter
    otherwise  
        err = MException('OPTI:SetFieldError','Unrecognized parameter name ''%s''.', field);
end
if(~isempty(err)), throw(err); end


function printfields(mode)
%Print out fields with defaults
if(~nargin), mode = 'ipopt'; end

fprintf('  IPOPT TERMINATION SETTINGS:\n');
fprintf('                       dual_inf_tol: [ Dual Infeasibility Tolerance {1} ]\n');
fprintf('                    constr_viol_tol: [ Absolute Constraint Violation Tolerance {1e-4} ]\n');
fprintf('                      compl_inf_tol: [ Absolute Complementarity Tolerance {1e-4} ]\n');
fprintf('                     acceptable_tol: [ Relative Acceptable Convergence Tolerance {1e-6} ]\n');
fprintf('                    acceptable_iter: [ Number of Acceptable Iterates before Termination {15} ]\n');
fprintf('         acceptable_constr_viol_tol: [ Absolute Acceptable Constraint Violation Tolerance {0.01} ]\n');
fprintf('            acceptable_dual_inf_tol: [ Acceptance Threshold for Dual Infeasibility {1e10} ]\n');
fprintf('           acceptable_compl_inf_tol: [ Absolute Acceptable Complementarity Tolerance {0.01} ]\n');
fprintf('          acceptable_obj_change_tol: [ Acceptable Objective Function Change Stopping Cirterion {1e20} ]\n');
fprintf('             diverging_iterates_tol: [ Threshold for Maximal Value of Primal Iterates {1e20} ]\n');

fprintf('\n  IPOPT NLP SCALING SETTINGS:\n');
fprintf('                 obj_scaling_factor: [ Scaling Factor for the Objective Function {1} ]\n');
fprintf('                 nlp_scaling_method: [ Scaling Technique {''gradient-based''}, ''none'' ]\n');
fprintf('           nlp_scaling_max_gradient: [ Maximum Gradient after NLP Scaling {100} ]\n');
fprintf('              nlp_scaling_min_value: [ Minimum Value of Gradient-Based Scaling Values {1e-8} ]\n');

fprintf('\n  IPOPT NLP SETTINGS:\n');
fprintf('                 bound_relax_factor: [ Factor for Initial Relaxation of Bounds {1e-8} ]\n');
fprintf('              honor_original_bounds: [ Indicates Whether Final Points should be Projected into Original Bounds ''no'', {''yes''} ]\n');
fprintf('       check_derivatives_for_naninf: [ Check Derivative Matrices for NaN/Inf {''no''}, ''yes'' ]\n');
fprintf('                nlp_lower_bound_inf: [ Bounds Less Than This Considered -Inf {-1e19} ]\n');
fprintf('                nlp_upper_bound_inf: [ Bounds Greater Than This Considered Inf {1e19} ]\n');
fprintf('           fixed_variable_treatment: [ Determines how Fixed Variables are Handled {''make_parameter''}, ''make_constraint'', ''relax_bounds'' ]\n');
fprintf('                     jac_c_constant: [ Indicates All Linear Equality Constraints {''no''}, ''yes'']\n');
fprintf('                     jac_d_constant: [ Indicates All Linear Inequality Constraints {''no''}, ''yes'' ]\n');
fprintf('                   hessian_constant: [ Indicates Quadratic Problem {''no''}, ''yes'' ]\n');

fprintf('\n  IPOPT INITIALIZATION SETTINGS:\n');
fprintf('                         bound_frac: [ Desired Minimum Relative Distance from Initial Point to Bound: {0.01} ]\n');
fprintf('                         bound_push: [ Desired Minimum Absolute Distance from Initial Point to Bound: {0.01} ]\n');
fprintf('                   slack_bound_frac: [ Desired Minimum Relative Distance from Initial Slack to Bound: {0.01} ]\n');
fprintf('                   slack_bound_push: [ Desired Minimum Absolute Distance from Initial Slack to Bound: {0.01} ]\n');
fprintf('                bound_mult_init_val: [ Initial Valye for All Bound Multipliers: {1} ]\n');
fprintf('               constr_mult_init_max: [ Maximum Allowed Least-Square Guess of Constraint Multipliers: {1000} ]\n');
fprintf('             bound_mult_init_method: [ Initialization Method for Bound Multipliers {''constant''}, ''mu-based'' ]\n');

fprintf('\n  IPOPT BARRIER PARAMETER SETTINGS:\n');
fprintf('                 mehrotra_algorithm: [ Use Mehrotra''s Algorithm (LP or Convex QP) {''no''}, ''yes'' ]\n');
fprintf('                        mu_strategy: [ Update Strategy for Barrier Parameter {''adaptive''}, ''monotone'' ]\n');
if(strcmpi(mode,'bonmin'))
    fprintf('                          mu_oracle: [ Oracle for New Parameter in Adaptive Strategy ''quality-function'', {''probing''}, ''loqo'' ]\n');
else
    fprintf('                          mu_oracle: [ Oracle for New Parameter in Adaptive Strategy {''quality-function''}, ''probing'', ''loqo'' ]\n');
end         
fprintf(' quality_function_max_section_steps: [ Maximum Number of Search Steps during Direct Search {8} ]\n');
fprintf('                    fixed_mu_oracle: [ Oracle for the Barrier Parameter When Fixed Mode ''quality_function'', ''probing'', ''loqo'', {''average_compl''} ]\n');
fprintf('          adaptive_mu_globalization: [ Globalization Strategy for the Adpative mu Selection Node {''obj-constr-filter''}, ''kkt-error'', ''never-monotone-mode'' ]\n');
fprintf('                            mu_init: [ Initial Value for the Barrier Parameter {0.1} ]\n');
fprintf('                        mu_max_fact: [ Factor for Initialization of Maximum Value for Barrier Parameter: {1000} ]\n');
fprintf('                             mu_max: [ Maximum Value for Barrier Parameter {100000} ]\n');
fprintf('                             mu_min: [ Minimum Value for Barrier Parameter {1e-11} ]\n');
fprintf('                          mu_target: [ Desired Value of Complementarity {0} ]\n');
fprintf('                 barrier_tol_factor: [ Factor for mu in Barrier Stop Test {10} ]\n');
fprintf('          mu_linear_decrease_factor: [ Determines Linear Decrease rate of Barrier Parameter {0.2} ]\n');
fprintf('      mu_superlinear_decrease_power: [ Determines Superlinear Decrease Rate of Barrier Parameter {1.5} ]\n');
        
fprintf('\n  IPOPT MULTIPLIER UPDATE SETTINGS:\n');
fprintf(['                        alpha_for_y: [ Method to Determine Step-Size for Constraint Multipliers {''primal''}, ''bound-mult'', ',...
        '''min'', ''max'', ''full'', ''min-dual-infeas'', ''safer-min-dual-infeas'', ''primal-and-full'', ''dual-and-full'', ''acceptor'' ]\n']);
fprintf('                    alpha_for_y_tol: [ Tolerance for Switching to Full Equality Multiplier Steps {10} ]\n');
fprintf('                           recalc_y: [ Tells Algorithm to Recalculate Constraint Multipliers as Least Square Estimates {''no''}, ''yes'' ]\n');  
fprintf('                  recalc_y_feas_tol: [ Feasibility Threshold for Recomputation of Multipliers {1e-6} ]\n');  

fprintf('\n  IPOPT LINE SEARCH SETTINGS:\n');
fprintf('                            max_soc: [ Maximum Number of Second Order Correction Steps {4} ]\n');
fprintf('    watchdog_shortened_iter_trigger: [ Number of Shortened Iterations that Trigger the Watchdog {10} ]\n');
fprintf('            watchdog_trial_iter_max: [ Maximum Number of Watchdog Iterations {3} ]\n');
fprintf('            accept_every_trial_step: [ Always Accept the First Trial Step {''no''}, ''yes'' ]\n');

fprintf('\n  IPOPT WARM START SETTINGS:\n');
fprintf('              warm_start_init_point: [ Warm-start for initial point {''no''}, ''yes'' ]\n');
fprintf('              warm_start_bound_push: [ Same as Bound Push for Regular Initializer {0.001} ]\n');
fprintf('              warm_start_bound_frac: [ Same as Bound Frac for Regular Initializer {0.001} ]\n');
fprintf('        warm_start_slack_bound_push: [ Same as Slack Bound Push for Regular Initializer {0.001} ]\n');
fprintf('        warm_start_slack_bound_frac: [ Same as Slack Bound Frac for Regular Initializer {0.001} ]\n');
fprintf('         warm_start_mult_bound_push: [ Same as Mult Bound Push for Regular Initializer {0.001} ]\n');
fprintf('           warm_start_mult_init_max: [ Maximum Initial Value for Equality Multipliers {1e6} ]\n');

fprintf('\n  IPOPT RESTORATION PHASE SETTINGS:\n');
if(strcmpi(mode,'bonmin'))
    fprintf('          expect_infeasible_problem: [ Enable Heuristics to Quickly Detect an Infeasible Problem ''no'', {''yes''} ]\n');
else
    fprintf('          expect_infeasible_problem: [ Enable Heuristics to Quickly Detect an Infeasible Problem {''no''}, ''yes'' ]\n');
end
fprintf('     expect_infeasible_problem_ctol: [ Threshold for Disabling "Expect Infeasible Problem" {0.001} ]\n');
fprintf('     expect_infeasible_problem_ytol: [ Multiplier Threshold for Activating "Expect Infeasible Problem" {1e8} ]\n');
fprintf('                   start_with_resto: [ Switch To Restoration Phase in First Iteration {''no''}, ''yes'' ]\n');
fprintf('soft_resto_pderror_reduction_factor: [ Required Reduction in Primal-Dual Error in the Soft Restoration Phase {0.9999} ]\n');
if(strcmpi(mode,'bonmin'))
    fprintf('   required_infeasibility_reduction: [ Required Reduction of Infeasibility before Leaving Restoration Phase {0.1} ]\n');
else
    fprintf('   required_infeasibility_reduction: [ Required Reduction of Infeasibility before Leaving Restoration Phase {0.9} ]\n');
end
fprintf('         bound_mult_reset_threshold: [ Threshold for Resetting Bound Multipliers after Restoration Phase {1000} ]\n');
fprintf('        constr_mult_reset_threshold: [ Threshold for Resetting Constraint Multipliers after Restoration Phase {0} ]\n');
fprintf('   evaluate_orig_obj_at_resto_trial: [ Determines if the Original Objective Function should be Evaluated at Restoration Trial Points ''no'', {''yes''} ]\n');

fprintf('\n  IPOPT LINEAR SOLVER SETTINGS:\n');
fprintf('                      linear_solver: [ Internal Linear System Solver {''MA57''}, ''PARDSIO'', ''MUMPS'' ]\n');
fprintf('              linear_system_scaling: [ Method for Scaling Linear System {''none''}, ''slack-based'' ]\n');
fprintf('           linear_scaling_on_demand: [ Indicate if Linear Scaling is Done On Demand (vs every linear system) ''no'', {''yes''} ]\n');
fprintf('               max_refinement_steps: [ Maximum Iterative Refinement Steps per Linear System Solve {10} ]\n');
fprintf('               min_refinement_steps: [ Minimum Iterative Refinement Steps per Linear System Solve {1} ]\n');

fprintf('\n  IPOPT HESSIAN PERTURBATION SETTINGS:\n');
fprintf('           max_hessian_perturbation: [ Maximum Regularization Parameter for Negative Curvature {1e20} ]\n');
fprintf('           min_hessian_perturbation: [ Smallest Perturbation of the Hessian Block {1e-20} ]\n');
fprintf('         first_hessian_perturbation: [ Size of first x-s Perturbation Tried {0.0001} ]\n');
fprintf('             perturb_inc_fact_first: [ Increase Factor for First x-s Perturbation {100} ]\n');
fprintf('                   perturb_inc_fact: [ Increase Factor for x-s Perturbation {8} ]\n');
fprintf('                   perturb_dec_fact: [ Decrease Factor for x-s Perturbation {0.333333} ]\n');
fprintf('      jacobian_regularization_value: [ Size of the Regularization for Rank-Deficient Constraint Jacobians {1e-8} ]\n');

fprintf('\n  IPOPT QUASI-NEWTON SETTINGS:\n');
fprintf('              hessian_approximation: [ Indicates Hessian Information to Use {''exact''}, ''limited-memory'' ]\n');
fprintf('         limited_memory_max_history: [ Maximum History for LBFGS Updates {6} ]\n');
fprintf('        limited_memory_max_skipping: [ Maximum Successive Iterations where Update is Skipped {2} ]\n');
fprintf('      limited_memory_initialization: [ Initialization Strategy for Quasi-Newton Update {''scalar1''}, ''scalar2'', ''scalar3'', ''scalar4'', ''constant'' ]\n');
fprintf('            limited_memory_init_val: [ Value for B0 in Low-Rank Update {1} ]\n');
fprintf('        limited_memory_init_val_max: [ Upper Bound on Value for B0 in Low-Rank Update {1e8} ]\n');
fprintf('        limited_memory_init_val_min: [ Lower Bound on Value for B0 in Low-Rank Update {1e-8} ]\n');
fprintf('   limited_memory_special_for_resto: [ Special Quasi-Newton for Restoration Phase {''no''}, ''yes'' ]\n');

fprintf('\n  IPOPT DERIVATIVE TEST SETTINGS:\n');
fprintf('                    derivative_test: [ Enable Derivative Checker {''none''}, ''first-order'', ''second-order'', ''only-second-order'' ]\n');
fprintf('       derivative_test_perturbation: [ Size of Finite Difference Perturbation {1e-8} ]\n');
fprintf('                derivative_test_tol: [ Threshold for Indicating Wrong Derivative {0.0001} ]\n');
fprintf('          derivative_test_print_all: [ Print all Derivative Test Information {''no''}, ''yes'' ]\n');
fprintf('        derivative_test_first_index: [ Index of First Quantity to be Checked by Derivative Checker {-2} ]\n');
fprintf('          point_perturbation_radius: [ Maximal Perturbation of an Evaluation Point {10} ]\n');

fprintf('\n  IPOPT MA57 LINEAR SOLVER SETTINGS:\n');
fprintf('                        ma57_pivtol: [ Pivot Tolerance {1e-8} ]\n');
fprintf('                     ma57_pivtolmax: [ Maximum Pivot Tolerance {0.0001} ]\n');
fprintf('                     ma57_pre_alloc: [ Safety Factor for Work Space Memory Allocation {1.05} ]\n');
fprintf('                   ma57_pivot_order: [ Controls Pivot Order {5} ]\n');
fprintf('             ma57_automatic_scaling: [ Controls Automatic Scaling ''no'', {''yes''} ]\n');
fprintf('                    ma57_block_size: [ Block size used by L3 BLAS {16} ]\n');
fprintf('             ma57_node_amalgamation: [ Node Amalgamation Parameter {16} ]\n');
fprintf('              ma57_small_pivot_flag: [ Small Entries are Removed and Correspondings Pivots Added {0} ]\n');

fprintf('\n  IPOPT MUMPS LINEAR SOLVER SETTINGS:\n');
fprintf('                       mumps_pivtol: [ Pivot Tolerance {1e-6} ]\n');
fprintf('                    mumps_pivtolmax: [ Maximum Pivot Tolerance {0.1} ]\n');
fprintf('                  mumps_mem_percent: [ Percent Increse in Estaimted Working Space {1000} ]\n');
fprintf('            mumps_permuting_scaling: [ Controls Permuting and Scaling {7} ]\n');
fprintf('                  mumps_pivot_order: [ Controls Pivot Order {7} ]\n');
fprintf('                      mumps_scaling: [ Controls Scaling {77} ]\n');

fprintf('\n  IPOPT PARDISO LINEAR SOLVER SETTINGS:\n');
fprintf(' pardiso_redo_symbolic_fact_only_if_inertia_wrong: [ Always redo symbolic factorization when elements were perturbed {''no''}, Only if Perturbed and Inertia Wrong (''yes'') ]\n');
fprintf('                      pardiso_order: [ Controls Fill-In Reduction Ordering Algorithm {''metis''}, ''amd'', ''pmetis'' ]\n');


fprintf('\n');
