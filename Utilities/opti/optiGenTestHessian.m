function testH = optiGenTestHessian(prob,x0)
%Generate a test hessian with dummy arguments

if(nargin < 2 || isempty(x0))
    % Get x0 (or try and make one)
    if(isfield(prob,'x0') && ~isempty(prob.x0))
        x0 = prob.x0;
    elseif(prob.sizes.ndec)
        x0 = opti_guessX0(prob,prob.sizes.ndec);
    elseif(~isempty(prob.A))
        x0 = opti_guessX0(prob,size(prob.A,2));
    elseif(~isempty(prob.nljacstr))
        ndec = size(prob.nljacstr(),2);
        x0 = opti_guessX0(prob,ndec);
    elseif(~isempty(prob.Hstr))
        ndec = size(prob.Hstr(),2);
        x0 = opti_guessX0(prob,ndec);
    else
        str = sprintf(['OPTI cannot determine the number of decision variables based on the arguments supplied.\n'...
               'Based on this, error checking cannot be performed, and the function cannot continue.\n\nPlease provide x0 to optiprob / opti in order to solve this problem.']);
        error(str); %#ok<SPERR> it doesn't take \n!
    end 
end

switch(nargin(prob.H))
    case 1
        try
            testH = prob.H(x0);
        catch ME
            fprintf(2,'There was an error when running a test Hessian function call. Please ensure this function exists and runs without error.\n\n');
            rethrow(ME);
        end        
    case 2
        try
            if(~isempty(prob.cl))
                v0 = ones(size(prob.cl));
            else
                v0 = ones(prob.sizes.nnlineq+prob.sizes.nnleq,1); %take a guess
            end
            testH = prob.H(x0,v0);
        catch ME
            fprintf(2,'There was an error when running a test Hessian function call. Please ensure this function exists and runs without error.\n\n');
            rethrow(ME);
        end        
    case 3
        %If loaded from an AMPL model, linear constraints also require a lambda entry (although make no difference to the Hessian)
        if(isfield(prob,'ampl') && ~isempty(prob.ampl.path))
            nlin = length(prob.rl);
        else
            nlin = 0;
        end
        try
            if(~isempty(prob.cl))
                v0 = ones(length(prob.cl)+nlin,1);
            else
                v0 = ones(prob.sizes.nnlineq+prob.sizes.nnleq+nlin,1); %take a guess
            end
            testH = prob.H(x0,1,v0);
        catch ME
            fprintf(2,'There was an error when running a test Hessian function call. Please ensure this function exists and runs without error.\n\n');
            rethrow(ME);
        end        
end