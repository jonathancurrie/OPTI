classdef opti < handle
%OPTI  Create an OPTI object for Linear, Quadratic & Nonlinear Optimization
%
%   optObj = opti('param1',value1,'param2',value2,...) creates an OPTI
%   object with the parameters 'param' set to their corresponding
%   values in 'value'. Parameters not specified will be set to the OPTI
%   default.
%
%   optObj = opti(...,'options',opts) creates the object with optiset
%   options supplied to the OPTI constructor.
%
%   optObj = opti(prob) creates an opti object based on the optimization
%   problem specified in prob. Use 'optiprob' to generate the structure
%   'prob'.
%
%   optObj = opti(prob,opts) creates the opti object with specified 
%   options. Use 'optiset' to generate the options structure 'opts'.
%
%   opti() prints a list of all possible fields and their function.
%
%
%   Typical usage examples:
%
%   1) LP
%   optObj = opti('f',f,'ineq',A,b,'bounds',lb,ub)
%
%   2) Bounded NLP
%   optObj = opti('obj',obj,'bounds',lb,ub,'x0',x0)
%
%   3) Option Setting
%   optObj = opti('obj',obj,'options',optiset('solver','nomad'))
%
%
%   For detailed documentation see the Wiki:
%      https://inverseproblem.co.nz/OPTI
%
%
%   See also opti.solve optiset opti.plot 
%
%   Copyright (C) 2011-2018 Jonathan Currie (www.inverseproblem.co.nz)
    
    properties (SetAccess = private)
        prob        % Problem Structure (optiprob)
        nlprob      % Nonlinear Problem Structure 
        opts        % Setup & Runtime options (optiset)
        sol         % Previous solution (Vector)
        obj         % Previous objective value
        ef          % Previous Exit Flag
        info        % Previous Info Structure
        fitstat     % Previous Fit Statistics Structure
    end
    
    methods
        
        %-- Constructor --%
        function op = opti(varargin) 
            %OPTI Constructor
            if(~nargin)
                optiprob;
                fprintf('\n OPTION SETTING:\n');
                fprintf('         options: [ Problem Setup & Solving Options specified via optiset() ]\n');
                fprintf('          solver: [ Optimization Solver ]\n'); 
                return;
            end
            [p,o] = opti.buildOpti(varargin{:});
            [op.prob,op.opts,op.nlprob] = opti.buildConfig(p,o);
            op.sol = []; op.ef = []; 
        end
        
        %-- Solve --%
        function [x,fval,exitflag,info] = solve(optObj,x0)
            %SOLVE  Solve OPTI Optimization Problem
            %
            %   [x,fval,exitflag,info] = solve(optObj) solves the optimization 
            %   problem specified by optObj and returns the solution vector x, 
            %   the function value at x in fval, an exitflag and information in 
            %   info regarding the solution.
            %
            %   [x,fval,exitflag,info] = solve(optObj,x0) solves the optimization
            %   problem using x0 as the initial guess. If x0 is specified within 
            %   opti() / optiprob(), this x0 will take precedence.
            %
            %   ExitFlags:
            %        1 - Converged / Terminated Successfully
            %        0 - Maximum Iterations / Function Evaluations / Time Exceeded
            %       -1 - Infeasible / Could Not Converge
            %       -2 - Unbounded / Solver Error
            %       -3 - Solver Specific Errors (Set option 'display' as 'iter')
            %       -5 - User Exited via Ctrl-C
            %
            %   Information Structure:
            %         Iterations - The number of iterations taken by the solver
            %               Time - The execution time of the solver as measured by MATLAB (tic + toc)
            %          Algorithm - The solver and algorithm being used
            %             Status - A status string indicating the solver specific exit message
            %             Lambda - A structure of the Lagrange multipliers at the solution (if available)
            
            %Check class isn't empty
            if(isempty(optObj.prob))
                error('You cannot solve an empty OPTI object!');
            end
            
            %Get Optional x0
            if(nargin < 2 || isempty(x0))
                if(~isempty(optObj.prob.x0) && ~any(isnan(optObj.prob.x0)))
                    x0 = optObj.prob.x0;
                else
                    x0 = [];
                end
            else
                optObj.prob.x0 = x0; %save into class
            end 
            %Solve Problem
            [x,fval,exitflag,info] = solveOpti(optObj,x0);
            %Save solve info
            optObj.sol = x;
            optObj.obj = fval;
            optObj.ef = exitflag;
            optObj.info = info;
            
            %Update ndec if empty
            if(isempty(optObj.prob.sizes.ndec)), optObj.prob.sizes.ndec = length(x0); end
            %Check for maximization
            if(optObj.prob.sense==-1), fval = -fval; end
            
            %Check if we should write an ampl solution file
            if(optObj.prob.ampl.writesol && asl('isopen') == 1)
                asl('writesol',info.Status(1:min(end,100)),x);
                %If the problem is not nonlinear, close the interface
                if(isempty(optObj.prob.fun) && isempty(optObj.prob.nlcon))
                    asl('close');
                end
            end
            
            %If less that 3 output arguments, user may not know if not
            %solved correctly. Depending on warning level, show a warning
            %to the user
            if(nargout < 3 && ~isempty(exitflag) && exitflag ~= 1 && optiWarnLevel(optObj.opts.warnings) >= 1 && isfield(info,'Status') && strcmpi(optObj.opts.display,'off'))
                optiwarn(['Your problem may not have solved correctly. Solver status:\n>> %s\n\nPlease enable solver output by changing the optiset setting ''display'' to ''iter'' to see more information.\n',...
                          'Alternatively you can avoid this warning by returning all four return arguments from solve(opti_obj).\n'],info.Status);
            end
        end
        
        
        %-- Multi-Solve --%
        function [x,fval,exitflag,info] = multisolve(optObj,x0,divisions,penalty,solveAll)
            %MULTISOLVE  Solve OPTI Global Optimization Problem with Multi-Start (beta)
            %
            %   x = multisolve(optObj) performs a search of the
            %   optimization problem before attempting to solve the problem
            %   at the best solution guesses found. The best optimized
            %   solution will be returned in x.
            %
            %   x = multisolve(optObj,x0) also includes a more detailed
            %   search around the initial guess x0.
            %
            %   x = multisolve(optObj,x0,divisions) determines the number
            %   of divisions each variable is searched within. The default
            %   will attempt to search around 100 points. For ND problems,
            %   a division of 10 will result in 10^N search points.
            %
            %   x = multisolve(optObj,x0,divisions,penalty) specifies the
            %   penalty factor applied to constraint violations. The
            %   default is 1e4. Larger values will be required for
            %   objective functions that exceed this value.
            %
            %   x = multisolve(optObj,x0,divisions,penalty,solveAll) 
            %   specifies at each search point an optimization is also 
            %   performed. Significantly longer to execute, but ensures the 
            %   search space is well covered!

            %Check class isn't empty
            if(isempty(optObj.prob))
                error('You cannot solve an empty OPTI object!');
            end
            
            %Optional args
            if(nargin < 5), solveAll = false; end
            if(nargin < 4), penalty = []; end
            if(nargin < 3), divisions = []; end
            
            %Get Optional x0
            if(nargin < 2 || isempty(x0))
                if(~isempty(optObj.prob.x0) && ~any(isnan(optObj.prob.x0)))
                    x0 = optObj.prob.x0;
                else
                    x0 = [];
                end
            else
                optObj.prob.x0 = x0; %save into class
            end
            
            %Solve Problem
            [x,fval,exitflag,info] = multiSolveOpti(optObj,x0,divisions,penalty,solveAll);
            %Save solve info into Object
            optObj.sol = x;
            optObj.obj = fval;
            optObj.ef = exitflag;
            optObj.info = info;
            
            %Update ndec if empty
            if(isempty(optObj.prob.sizes.ndec)), optObj.prob.sizes.ndec = length(x0); end
            %Check for maximization
            if(optObj.prob.sense==-1), fval = -fval; end
            
            %Check if we should write an ampl solution file
            if(optObj.prob.ampl.writesol && asl('isopen') == 1)
                asl('writesol',info.Status(1:min(end,100)),x);
                %If the problem is not nonlinear, close the interface
                if(isempty(optObj.prob.fun) && isempty(optObj.prob.nlcon))
                    asl('close');
                end
            end
        end        
        
        %-- Check Solution --%
        function [ok,msg] = checkSol(optObj,tol)
            %CHECKSOL  Check Optimization Solution for Errors
            %
            %   [ok,msg] = checkSol(optObj) checks the solution stored in
            %   optObj for error messages or broken constraints. It returns
            %   true in 'ok' if no problems detected, otherwise it returns
            %   false with the problem in 'msg'.
            %
            %   [ok,msg] = checkSol(optObj,tol) checks with respect to a
            %   given tolerance. The default is 1e-6.
            
            if(isempty(optObj.sol))
                error('This OPTI object has not been solved yet!');
            end
            if(nargin < 2), tol = 1e-6; end
            
            [ok,msg] = checkOptiSol(optObj,tol);
        end  
        
        %-- Display --%
        function display(optObj)               
            %Check if empty
            if(isempty(optObj.prob))
                disp(' ');
                disp('------------------------------------------------------');
                disp('Empty OPTI Object');
                disp('------------------------------------------------------');
            else
                %Otherwise normal method
                displayOPTI(optObj);
            end
        end
            
        %-- Plot --%
        function plot(optObj,scale,doLog,npts)
            %PLOT  Plot optimization problem (1D-5D only)
            %
            %   plot(optObj) plots the optimization field for the current
            %   OPTI object.
            %
            %   plot(optObj,scale) plots with a user defined zoom level.
            %   scale is defined as the range +- of the solution to be
            %   drawn OR if supplied as a vector controls the bounds on the
            %   plot [x1min x1max x2min x2max ... xNmin xNmax]
            %
            %   plot(obtObj,scale,dolog) plots the log of the objective
            %   function (NL only)
            %
            %   plot(obtObj,scale,dolog,npts) controls the number of points
            %   used for generating the objective and constraint contours
            
            %Check class isn't empty
            if(isempty(optObj.prob))
                error('You cannot plot an empty OPTI object!');
            end
                        
            %Get Optional Arguments
            if(nargin < 4), npts = []; end
            if(nargin < 3), doLog = false; end
            if(nargin < 2), scale = []; end

            %Check 2D for problems other than least squares fits
            if(~any(strcmpi(optObj.prob.type,{'NLS','DNLS'})) && optObj.prob.sizes.ndec > 5)
                error('You can only plot problems with one - five decision variables!');
            end
            
            %Plot
            plotOptiProb(optObj.prob,optObj.opts,optObj.sol,optObj.fitstat,scale,doLog,npts,'normal');                    
        end      
        
        
        %-- Multi-Plot--%
        function multiplot(optObj,doLog,npts)
            %MULTIPLOT  Plot optimization problem and overlay multi-start search area
            %
            %   multiplot(optObj) plots the optimization field for the current
            %   OPTI object and overlays the multi-start search area
            %
            %   multiplot(obtObj,dolog) plots the log of the objective
            %   function (NL only)
            %
            %   multiplot(obtObj,dolog,npts) controls the number of points
            %   used for generating the objective and constraint contours
            
            if(nargin < 3), npts = []; end
            if(nargin < 2), doLog = false; end
            
            %Check class isn't empty
            if(isempty(optObj.prob))
                error('You cannot plot an empty OPTI object!');
            end
            %Check 2D
            if(optObj.prob.sizes.ndec > 5)
                error('You can only multiplot problems with one - five decision variables!');
            end
            
            %Plot         
            plotOptiProb(optObj.prob,optObj.opts,optObj.sol,optObj.fitstat,[],doLog,npts,'multi'); 
        end
        
        
        
        
        %-- Stats --%
        function [StatsStruct] = fitStats(optObj,limit,printStats)
            %CONFIDENCE  Calculate approximate ANOVA and Confidence Statistics for NLS/DNLS Fit
            %
            %   StatsStruct = fitStats(optObj) returns a structure of
            %   statistics about the fit, as detailed below:
            %       SSE - Sum Of Squared Errors
            %       Rsquare - R^2 Value
            %       AdjRSquare - Adjusted R^2 Value
            %       RMSE - Root Mean Squared Error
            %       DFE - Degrees of Freedom
            %       ContInt - Confidence Interval (theta +- ConfInt)
            %       ConfBnds - Confidence Bounds (ypred +- ConfBnds)
            %       Cov - System Covariance Matrix
            %       BndIdx - Logical Array where 1s indicate active bounds
            %       at the solution (and therefore dropped from Cov)
            %       Param - Parameter Statistics Structure with fields:
            %           StdError - Standard Error
            %           tStat - t Statistic
            %           pValues - p Values
            %       Model - Model Statistics Stucture with fields:
            %           Fstat - F Statistic
            %           pValue - p value
            %           Int - Whether Model has an Intercept or Not
            %       Conf - User specified Confidence limit                   
            %
            %   For a description of each statistic, consult the following
            %   site:
            %       http://www.ats.ucla.edu/stat/sas/output/reg.htm
            %
            %   StatsStruct = fitStats(optObj,limit) allows the user to 
            %   specify the confidence limit (0.0 < limit < 1.0), 
            %   e.g. 0.95 = 95%
            %
            %   fitStats(optObj,limit,printStats) optionally prints the
            %   statistics to the command window. Defaults if no output is
            %   selected.
            
            %Check class isn't empty or unsolved
            if(isempty(optObj.prob) && isempty(optObj.sol))
                error('You cannot calculate statistics of an empty OPTI object!');
            end
            if(isempty(optObj.sol))
                error('You cannot calculate statistics of an unsolved OPTI object!');
            end
                        
            %Get Optional Arguments
            if(nargin < 3), printStats = false; end
            if(nargin < 2), limit = []; end            
            
            %Calculate
            StatsStruct = calcStatistics(optObj,limit);   
            %Save to Object
            optObj.fitstat = StatsStruct;
%             %Optionally Plot Ellipse
%             if(plotEllipse)
%                 if(length(optObj.sol) < 2 || length(optObj.sol) > 3)
%                     optiwarn('opti:confEllipse','Cannot plot confidence ellipse/ellipsoid - only 2/3 parameter problems are supported');
%                 else
%                     plotConfEllipse(optObj.sol,StatsStruct);
%                 end
%             end
            %If no output arguments, display fit statistics
            if(nargout < 1 || printStats)
               	optiPrintFitStats(optObj.fitstat,optObj.sol);
            end
        end 
        
        
        %-- write --%
        function write(optObj,filename,type)
            %WRITE  Write an OPTI Problem to a MPS/QPS/LP/SDPA/SEDUMI/GAMS file
            %
            %   write(optObj,filename) writes the OPTI problem to a
            %   filetype specified by the extension of filename.
            %
            %   write(optObj,filename,type) writes the OPTI problem to a
            %   filetype specified by the type argument.
            %
            %   Valid File Types & Extensions:
            %       Linear & Quadratic Problems
            %           - MPS [.mps] (LP/MILP)
            %           - QPS [.qps] (QP/MIQP)
            %           - LP  [.lp]  (LP)
            %       Semidefinite Problems
            %           - SDPA, DAT [.dat] (Dense)
            %           - SDPA-S, DAT-S [.dat-s] (Sparse)
            %           - MAT [.mat] (SeDuMi Sparse)
            %       General
            %           - GMS [.gms] (GAMS)
            
            %Check class isn't empty
            if(isempty(optObj.prob)), error('You cannot write an empty OPTI object!'); end            
            if(nargin < 3), type = []; end
            %Check if we have a type
            if(~isempty(type))
                switch(lower(type))
                    case {'sdpa','sdpa-s','sedumi','sdpas','sdp','dat','dat-s','mat'}
                        sdpWrite(optObj.prob,filename,type);
                        return;
                    case {'mps','qps','lp'}            
                        coinWrite(optObj.prob,filename,type);
                        return;
                    case 'gms'
                        gamsWrite(optObj.prob,filename);
                    otherwise
                        error('Unknown file type to write to!');
                end
            end
            %If not, go on extension
            if(isempty(strfind(filename,'.')))
                error('You must supply a file extension if you do not specify a file type!');
            end
            ext = regexp(filename,'\.','split');
            switch(lower(ext{end}))
                case {'dat','mat','dat-s'}
                    sdpWrite(optObj.prob,filename,type);
                case {'mps','qps','lp'}            
                    coinWrite(optObj.prob,filename,type);
                case 'gms'            
                    gamsWrite(optObj.prob,filename);
                otherwise
                    error('Unknown file extension to write to! For non-standard file extensions please specify the file type.');
            end
        end
        
        
        %-- getProb --%
        function retprob = getProb(optObj)
            %GETPROB  Return optiprob compatible problem structure
            %
            %   prob = getProb(optObj)
            
            %Check class isn't empty
            if(isempty(optObj.prob))
                error('You cannot get the problem from an empty OPTI object!');
            end
            
            retprob = optObj.prob;
            retprob = rmfield(retprob,{'iscon','sizes','numdif','type','ampl','save'});
            retprob.int = retprob.int.str;
        end      
        
        %-- Set Solution (unit testing only use) --%
        function setSolution(optObj, fval, sol)
            %Check class isn't empty
            if(isempty(optObj.prob))
                error('You cannot get the problem from an empty OPTI object!');
            end
            optObj.obj = fval;
            if (numel(sol) ~= optObj.prob.sizes.ndec)
                error('Solution vector is the wrong size');
            else
                optObj.sol = sol(:);
            end
            optObj.info = struct('Status','HARDCODED SOLUTION','Algorithm','USER SPECIFIED');
        end
    end
    
    methods (Static)
        [prob,opts] = buildOpti(varargin);
        [prob,opts,nlprob] = buildConfig(prob,opts);
    end
    
end

