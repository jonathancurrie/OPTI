function options = nomadset(varargin)
%NOMADSET  Create or alter the options for Optimization with NOMAD
%
% options = nomadset('param1',value1,'param2',value2,...) creates an
% NOMAD options structure with the parameters 'param' set to their 
% corresponding values in 'value'. Parameters not specified will be set to 
% the NOMAD default.
%
% options = nomadset(oldopts,'param1',value1,...) creates a copy of the old
% options 'oldopts' and then fills in (or writes over) the parameters
% specified by 'param' and 'value'.
%
% options = nomadset() creates an options structure with all fields set to
% NOMADSET defaults.
%
% nomadset() prints a list of all possible fields and their function.

%   Copyright (C) 2012 Jonathan Currie (I2C2)

%Valid direction types
global dirtypes
dirtypes = {'ortho 1','ortho 2','ortho n+1','ortho n+1 quad','ortho n+1 neg','ortho','ortho 2n','lt 1','lt 2','lt 2n','lt n+1','gps binary','gps 2n static',...
            'gps 2n rand','gps n+1 static uniform','gps n+1 static','gps n+1 rand uniform','gps n+1 rand'};

% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    printfields();
    return
end

%Names and Defaults
Names = {'bb_input_type','bb_output_type','direction_type','f_target','initial_mesh_size','lh_search','max_bb_eval',...
         'max_time','model_eval_sort','model_search','multi_nb_mads_runs','multi_overall_bb_eval','opportunistic_eval','opportunistic_lh','seed','vns_search',...
         'cache_search','h_max_0','h_min','h_norm','initial_mesh_index','l_curve_target','max_cache_memory','max_consecutive_failed_iterations',...
         'max_eval','max_iterations','max_mesh_index','max_sim_bb_eval','mesh_coarsening_exponent','mesh_refining_exponent','mesh_update_basis',...
         'min_mesh_size','min_poll_size','model_eval_sort_cautious','model_search_max_trial_pts','model_search_optimistic','model_search_proj_to_mesh',...
         'model_quad_max_y_size','model_quad_min_y_size','model_quad_radius_factor','model_quad_use_wp','multi_f_bounds','multi_formulation','multi_use_delta_crit',...
         'opportunistic_cache_search','opportunistic_lucky_eval','opportunistic_min_eval','opportunistic_min_f_imprvmt','opportunistic_min_nb_success','rho','scaling','sec_poll_dir_type',...
         'snap_to_bounds','speculative_search','stat_sum_target','stop_if_feasible','add_seed_to_file_names','cache_file','cache_save_period','display_degree','display_all_eval',...
         'history_file','solution_file','stats_file','param_file','iterfun','disable','epsilon','opt_only_sgte','sgte_cost','sgte_eval_sort','has_sgte',...
         'sgte_cache_file','max_sgte_eval'}';
Defaults = {[],[],'ortho n+1 quad',[],[],[],[],[],1,1,[],[],1,1,0,0,0,1e20,[],'L2',[],[],2000,[],[],[],[],[],1,-1,4,... %mesh_update_basis
            [],[],0,10,1,1,500,[],2,0,[],'product',0,0,[],[],[],-1,0.1,[],[],1,1,[],0,1,[],25,0,0,[],[],...
            [],[],[],[],1e-13,0,[],1,0,[],[]}';        

%Collect Sizes and lowercase matches         
m = size(Names,1); numberargs = nargin;
%Create structure with all names and default values
st = [Names,Defaults]'; options = struct(st{:});

% Check we have char or structure input. If structure, insert arguments
i = 1;
while i <= numberargs
    arg = varargin{i};
    if ischar(arg)
        break;
    end
    if ~isempty(arg)
        if ~isa(arg,'struct')
            error('An argument was supplied that wasn''t a string or struct!');
        end
        for j = 1:m
            if any(strcmp(fieldnames(arg),Names{j,:}))
                val = arg.(Names{j,:});
            else
                val = [];
            end
            if(~isempty(val))
                checkfield(Names{j,:},val);
                options.(Names{j,:}) = val;
            end
        end
    end
    i = i + 1;
end

%Check we have even number of args left
if rem(numberargs-i+1,2) ~= 0
    error('You do not have a value specified for each field!');
end

%Go through each argument pair and assign to correct field
expectval = 0; %first arg is a name
while i <= numberargs
    arg = varargin{i};

    switch(expectval)
        case 0 %field
            if ~ischar(arg)
                error('Expected field name to be a string! - Argument %d',i);
            end
            j = find(strcmp(arg,Names) == 1);
            if isempty(j)  % if no matches
                error('Unrecognised parameter %s',arg);
            elseif(length(j) > 1)
                error('Ambiguous parameter %s',arg);
            end
            expectval = 1; %next arg is a value
        case 1
            if(~isempty(arg))
                if(iscell(arg))
                    checkfield(Names{j,:},arg);
                    options.(Names{j,:}) = arg;
                else
                    if(ischar(arg)), arg = lower(arg); end
                    checkfield(Names{j,:},arg);
                    options.(Names{j,:}) = arg;
                end
            end
            expectval = 0;
    end
    i = i + 1;
end

if expectval %fallen off end somehow
    error('Missing value for %s',arg);
end
%Ensure we have absolute path
if(~isempty(options.param_file))
    options = checkParamFile(options);
end


function checkfield(field,value)
%Check a field contains correct data type
global dirtypes
err = [];

switch lower(field)   
    %Scalar double
    case {'f_target','max_bb_eval','max_time','model_eval_sort',...
         'model_search','multi_nb_mads_runs','multi_overall_bb_eval','opportunistic_eval','opportunistic_lh','seed','vns_search',...
         'cache_search','h_max_0','h_min','initial_mesh_index','l_curve_target','max_cache_memory','max_consecutive_failed_iterations',...
         'max_eval','max_iterations','max_mesh_index','max_sim_bb_eval','mesh_coarsening_exponent','mesh_refining_exponent','mesh_update_basis',...
         'model_eval_sort_cautious','model_search_max_trial_pts','model_search_optimistic','model_search_proj_to_mesh',...
         'model_quad_max_y_size','model_quad_min_y_size','model_quad_radius_factor','model_quad_use_wp','multi_use_delta_crit',...
         'opportunistic_cache_search','opportunistic_lucky_eval','opportunistic_min_f_imprvmnt','opportunistic_min_nb_success','rho',...
         'snap_to_bounds','speculative_search','stat_sum_target','stop_if_feasible','add_seed_to_file_names','cache_save_period','display_degree',...
         'display_all_eval','has_sgte','sgte_cost','sgte_eval_sort','opt_only_sgte','epsilon','max_sgte_eval'}
        if(~isscalar(value) || ~isnumeric(value) || ~isreal(value) || ~isa(value,'double'))
            err = MException('NOMAD:SetFieldError','Parameter ''%s'' should be a real scalar double',field);
        end
    
    %Direction Type
    case {'direction_type','sec_poll_dir_type'}
        if(~ischar(value) || ~any(strcmp(value,dirtypes)))
            err = MException('NOMAD:SetFieldError','Parameter ''%s'' should match one of those displayed by nomadset()',field);
        end   
        
    %Disable
    case 'disable'
        if(~ischar(value) || ~any(strcmp(value,{'models'})))
            err = MException('NOMAD:SetFieldError','Parameter ''%s'' should be one of {''models''}',field);
        end 
    
    %String
    case {'cache_file','history_file','solution_file','stats_file','param_file','multi_formulation','h_norm',...
          'lh_search','multi_f_bounds','sgte_cache_file'}
        if(~ischar(value))
            err = MException('NOMAD:SetFieldError','Parameter ''%s'' should be a char array (string)',field);
        end 
     
    %string or cell of strings
    case {'bb_output_type','bb_input_type','initial_mesh_size','min_mesh_size','min_poll_size','scaling'}
        if(~ischar(value) && ~iscell(value))
            err = MException('NOMAD:SetFieldError','Parameter ''%s'' should be a string or cell array of strings',field);
        elseif(iscell(value))
            for i = 1:numel(value)
                if(~ischar(value{i}))
                    err = MException('NOMAD:SetFieldError','Parameter ''%s'' cell(%d) should be a string',field,i);
                    break;
                end
            end
        end
    
    %function_handle
    case 'iterfun'
        if(~isa(value,'function_handle'))
            err = MException('NOMAD:SetFieldError','Parameter ''%s'' should be a function handle',field);
        end
        
    otherwise  
        err = MException('OPTI:SetFieldError','Unrecognized parameter name ''%s''.', field);
end
if(~isempty(err)), throw(err); end


function opts = checkParamFile(opts)
%Check parameter file exists and collect absolute path
filename = opts.param_file;
%See if we got an absolute path
switch(computer)
    case {'PCWIN','PCWIN64'}
        if(~isempty(strfind(filename,':')))
            p = filename;
            %Check it exists
            if(exist(filename,'file') ~= 2)
                error('Cannot locate absolute file %s!',filename);
            end
        else
            p = which(filename);
            %Check it exists
            if(isempty(p))
                error('Cannot locate file %s!',filename);
            end
        end
    otherwise %Linux / Mac
        if(filename(1) == '/') %absolute path (I think)
            p = filename;
            %Check it exists
            if(exist(filename,'file') ~= 2)
                error('Cannot locate absolute file %s!',filename);
            end
        else
            p = which(filename);
            %Check it exists
            if(isempty(p))
                error('Cannot locate file %s!',filename);
            end
        end
end
%Save full path
opts.param_file = p;


function printfields()
%Print out fields with defaults

global dirtypes

fprintf('\nBASIC PROBLEM PARAMETERS:\n');
fprintf('              bb_input_type: [ Blackbox input types, specified as a cell array of strings, overwrites xtype {[]} ] \n');
fprintf('             bb_output_type: [ Blackbox output types, specified as a cell array of strings (see below) ] \n');

fprintf('\nBASIC ALGORITHMIC PARAMETERS:\n');
fprintf('             direction_type: [ Type of directions for the poll (see below) {''ortho n+1 quad''} ] \n');
fprintf('                   f_target: [ Stop if f(x) < ftarget {[]} ] \n');
fprintf('          initial_mesh_size: [ Initial Mesh Size (del_0^m) {[]} ]\n');
fprintf('                  lh_search: [ Latin-Hypercube Search (po and pi) {[]} ]\n');
fprintf('                max_bb_eval: [ Maximum number of blackbox evaluations {[]} ] \n');
fprintf('                   max_time: [ Maximum execution time {[]} ] \n');

fprintf('\nADVANCED ALGORITHMIC PARAMETERS:\n');
fprintf('                  sgte_cost: [ The cost of c surrogate evaluations is equivalent to one blackbox evaluation: {[]} ] \n');
fprintf('             sgte_eval_sort: [ If surrogates are used to sort list of trial points: Off (0), On {1} ] \n');
fprintf('               cache_search: [ Use Cache search: Off {0}, On (1) ] \n');
fprintf('                    disable: [ Forcefully disable NOMAD features {[]} ] \n');
fprintf('                    h_max_0: [ Initial value of hmax {1e20} ] \n');
fprintf('                      h_min: [ x is feasible if h(x) >= v {0.0} ] \n');
fprintf('                     h_norm: [ Norm used to compute h: ''L1'', {''L2''}, ''Linf'' ] \n');
fprintf('                   has_sgte: [ Indicates if the problem has a surrogate function: No {0}, Yes {1} ] \n');
fprintf('         initial_mesh_index: [ Initial Mesh Index {0} ] \n');
fprintf('             l_curve_target: [ NOMAD terminates if objective may not reach this value {[]} ] \n');
fprintf('           max_cache_memory: [ Maximum cache memory {2000} MB ] \n');
fprintf('   max_consecutive_failed_iterations: [ Maximum number of failed MADS iterations {[]} ] \n');
fprintf('                   max_eval: [ Maximum number of evaluations (includes cache and bb) {[]} ] \n');
fprintf('             max_iterations: [ Maximum number of iterations {[]} ] \n');
fprintf('             max_mesh_index: [ Maximum Mesh Index {[]} ] \n');
fprintf('              max_sgte_eval: [ Maximum Number of Surrogate Evaluations {[]} ] \n');
fprintf('            max_sim_bb_eval: [ Maximum Simulated BB evaluations {[]} ] \n');
fprintf('   mesh_coarsening_exponent: [ Mesh Coarsening Exponent {1} ] \n');
fprintf('     mesh_refining_exponent: [ Mesh Refining Exponent {-1} ] \n');
fprintf('          mesh_update_basis: [ Mesh Update Basis {4} ] \n');
fprintf('              min_mesh_size: [ Minimum Mesh Size (del_min^m) {[]} ] \n');
fprintf('              min_poll_size: [ Minimum Poll Size (del_min^p) {[]} ] \n');
fprintf('            model_eval_sort: [ Use Quadratic Model to order points: Off (0), On {1} ] \n');
fprintf('               model_search: [ Use Quadratic Model searches: Off (0), On {1} ] \n');
fprintf('    model_search_optimistic: [ If the model search is optimistic : No (0), Yes {1} ] \n');
fprintf('             multi_f_bounds: [ Vector of [f1min, f1max, f2min, f2max] for biobjective surf calculation {[]} ] \n');
fprintf('         multi_nb_mads_runs: [ Number of MADS runs in Biobjective Optimization {[]} ] \n');
fprintf('      multi_overall_bb_eval: [ Maximum number of BB evals for all MADS runs {[]} ] \n');
fprintf(' opportunistic_cache_search: [ Opportunistic cache search: Off {0}, On (1) ] \n');
fprintf('         opportunistic_eval: [ Opportunistic strategy: Off (0), On {1} ] \n');
fprintf('           opportunistic_lh: [ Opportunistic strategy for LH search: Off(0), On {1} ] \n');
fprintf('     opportunistic_min_eval: [ Do not terminate below i evaluations {[]} ] \n');
fprintf('                        rho: [ Parameter of the progressive barrier {0.1} ] \n');
fprintf('                    scaling: [ Scaling on the variables (vector of scaling values for each variable) {[]} ] \n');
fprintf('                       seed: [ Random Seed: {0} ] \n');
fprintf('             snap_to_bounds: [ Snap to boundary trial points that are generated outside of bounds: No (0), Yes {1} ] \n');
fprintf('         speculative_search: [ MADS speculative search: No (0), Yes {1} ] \n');
fprintf('            stat_sum_target: [ Terminates if stat_sum reaches this value {[]} ] \n');
fprintf('           stop_if_feasible: [ Stop on first feasible solution: Off {0}, On (1) ] \n');
fprintf('                 vns_search: [ Variable Neighbourhood Search (Multiple Minima Problems): Off {0.0}, Max Searches (1.0) ] \n');

fprintf('\nDEVELOPER PARAMETERS:\n');
fprintf('                    epsilon: [ Precision on real numbers {1e-13} ]\n');
fprintf('   model_eval_sort_cautious: [ If the model ordering strategy is cautious: No {0}, Yes (1) ] \n');
fprintf(' model_search_max_trial_pts: [ Maximum trial points for one model search {10} ] \n');
fprintf('  model_search_proj_to_mesh: [ If model search trial points are projected to mesh: No (0), Yes {1} ] \n');
fprintf('      model_quad_max_y_size: [ Upper limit on the size of interpolation sets for quadratic models {500} ] \n');
fprintf('      model_quad_min_y_size: [ Inf limit on the size of interpolation sets for quadratic models {[]} ] \n');
fprintf('   model_quad_radius_factor: [ Quadratic Model search radius factor {2.0} ] \n');
fprintf('          model_quad_use_wp: [ Enable strategy to maintain well-poised quadratic models: Off {0}, On (1) ] \n');
fprintf('          multi_formulation: [ Single objective reformulation: ''normalized'', {''product''}, ''dist_l1'', ''dist_l2'', ''dist_linf'' ] \n');
fprintf('       multi_use_delta_crit: [ Use stopping criterion based on the delta criterion: Off {0}, On (1) ] \n');
fprintf('   opportunistic_lucky_eval: [ Perform an additional BB eval after an improvement: Off {0}, On (1) ] \n');
fprintf('  opportunistic_min_f_imprvmt: [ Terminate only if f is reduced by r%% {[]} ] \n');
fprintf(' opportunistic_min_nb_success: [ Do not terminate before i successes {-1} ] \n');
fprintf('              opt_only_sgte: [ Minimize only with surrogates: No {0}, Yes (1) ]\n');
fprintf('          sec_poll_dir_type: [ Type of directions for the secondary poll (see below) {[]} ] \n');

fprintf('\nOUTPUT PARAMETERS:\n');
fprintf('     add_seed_to_file_names: [ If the seed is added to the output file names: Off (0), On {1} ]\n');
fprintf('                 cache_file: [ Cache file: Off {[]}, On ''filename'' ]\n');
fprintf('          cache_save_period: [ Cache files are saved every i iterations: {25} ]\n');
fprintf('             display_degree: [ Display level: Increasing (>=0 and <=3) {0}  ] \n');
fprintf('           display_all_eval: [ Display all Points: Off {0}, On (1) ] \n');
fprintf('               history_file: [ File containing all trial points: Off {[]}, On ''filename'' ]\n');
fprintf('              solution_file: [ File to save the current best feasible point: Off {[]}, On ''filename'' ]\n');
fprintf('                 stats_file: [ File for screen dump: Off {[]}, On ''filename'' ]\n');
fprintf('            sgte_cache_file: [ Surrogate cache file: Off {[]}, On ''filename'' ]\n');

fprintf('\nMEX INTERFACE PARAMETERS:\n');
fprintf('                 param_file: [ Read NOMAD options from a parameters file (overwrites these options): Off {[]}, On ''filename'' ]\n');
fprintf('                    iterfun: [ Callback function run every function evaluation: Off {[]}, On function_handle ]\n');


fprintf('\n\nValid Direction Types (string):\n');
for i = 1:length(dirtypes)
    fprintf('''%s''\n',dirtypes{i});
end
fprintf('\nValid BB_OUTPUT_TYPE types (cell array of strings):\n');
fprintf('OBJ: Objective\n');
fprintf(' PB: Progressive Barrier (default constraint type)\n');
fprintf(' EB: Extreme Barrier\n');
fprintf('PEB: Hybrid Constraint (PB/EB)\n');
fprintf('  F: Filter\n');




