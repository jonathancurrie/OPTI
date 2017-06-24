function ok = opti_Install_Test(verb)
% Installation Test File
%
% A collection of tests to verify the installation of the toolbox has been
% successful. 

%   Copyright (C) 2012 Jonathan Currie (IPL)

% 14/3/12 Changed to Rel Error Check
% 8/2/14 Changed tol to 1e-3

%Set Verbosity
if(nargin < 1)
    verb = 1; %default show all messages
end

%Set default ok
ok = 1;

fprintf('\nChecking OPTI Toolbox Installation:\n');

%% TEST 1 - Check Main Paths
if(verb); fprintf('Checking Paths...                '); end
paths = {'@opti','Utilities','Help','Solvers','opti'};
len = length(paths);
for i = 1:len
    pass = exist(paths{i},'dir');
    if(~pass)
        fprintf('\nFailed Path Check on %s\n',paths{i});
        ok = 0;
        return
    end
end
if(verb); fprintf('Ok\n'); end

%% Setup tests to run
tests = {'lp','milp','qp','miqp','sdp','nls','nlp','minlp'};
tres = ones(size(tests));

%Loop for each test set
for t = 1:length(tests)    
    load([tests{t} '_test_results']); %load pre solved results 
    probs = eval([tests{t} '_tprob']);
    sols = eval([tests{t} '_sval']);
    msolvers = optiSolver(tests{t}); 
    ind = strcmpi(msolvers,'MATLAB'); msolvers(ind) = [];    
    ind = strcmpi(msolvers,'GMATLAB'); msolvers(ind) = [];
    %Check and remove solvers for certain problems
    switch(lower(tests{t}))
        case 'nls'
            ind = strcmpi(msolvers,'NOMAD'); msolvers(ind) = [];
            ind = strcmpi(msolvers,'PSWARM'); msolvers(ind) = [];
            ind = strcmpi(msolvers,'NLOPT'); msolvers(ind) = []; 
        case 'nlp'
            ind = strcmpi(msolvers,'NLOPT'); msolvers(ind) = [];
            ind = strcmpi(msolvers,'PSWARM'); msolvers(ind) = [];
            ind = strcmpi(msolvers,'SCIP'); msolvers(ind) = []; %problems with trig functions
            ind = strcmpi(msolvers,'BARON'); msolvers(ind) = []; %problems with trig functions
        case 'lp'
            ind = strcmpi(msolvers,'DSDP'); msolvers(ind) = []; %no good with equalities
            ind = strcmpi(msolvers,'LIPSOL'); msolvers(ind) = []; %no good with infinite lower bounds
        case 'minlp'
            ind = strcmpi(msolvers,'BARON'); msolvers(ind) = []; %problems with trig functions
    end            
    noS = length(msolvers); 
    if(noS == 0)
        continue;
    end    
    noP = length(probs);
    if(verb); fprintf('Checking %5s Solver Results... ',upper(tests{t})); end
    for i = 1:noS
        for j = 1:noP
            if (strcmpi(msolvers{i},'nomad') && j == 1)
                continue; %annoying problems with v3.8.1 on this
            end
            Opt = opti(probs{j},optiset('solver',msolvers{i},'warnings','off'));
            [~,fval] = solve(Opt);
            if(sols(j) == 0)
                if(abs(fval-sols(j)) > 1e-3)
                    fprintf('\nFailed %s Result Check on Test %d, Solver: %s',tests{t},j,msolvers{i});
                    tres(t) = 0;
                end
            else
                if(abs((fval-sols(j))/fval) > 1e-3) 
                    fprintf('\nFailed %s Result Check on Test %d, Solver: %s',tests{t},j,msolvers{i});
                    tres(t) = 0;
                end
            end
        end
    end
    if(verb && tres(t)), fprintf('Ok\n'); else, fprintf('\n\n'); end 
end

%% Final
if(ok && all(tres))
    ok = 1;
    fprintf('\nToolbox Checked Out Ok! - Enjoy\n');
else
    ok = 0;
    fprintf(2,'\nThere was an error during post installation checks, please ensure your download was not corrupted!\n');
end

end

