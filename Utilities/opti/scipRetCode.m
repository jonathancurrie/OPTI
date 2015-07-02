function [status,exitflag] = scipRetCode(code)
%Returns Status associated with SCIP Return Code

%Note v3.2.0 some numbers changed! Below as per scip.m
%   Return Status:
%       0 - Unknown
%       1 - User Interrupted
%       2 - Node Limit Reached
%       3 - Total Node Limit Reached
%       4 - Stall Node Limit Reached
%       5 - Time Limit Reached
%       6 - Memory Limit Reached
%       7 - Gap Limit Reached
%       8 - Solution Limit Reached
%       9 - Solution Improvement Limit Reached
%      10 - Restart Limit Reached
%      11 - Problem Solved to Optimality
%      12 - Problem is Infeasible
%      13 - Problem is Unbounded
%      14 - Problem is Either Infeasible or Unbounded

switch(code)    
    case 1, status = 'User Exited'; exitflag = -5;
    case 2, status = 'Node Limit Reached'; exitflag = 0;
    case 3, status = 'Total Node Limit Reached'; exitflag = 0;
    case 4, status = 'Stall Node Limit Reached'; exitflag = 0;
    case 5, status = 'Time Limit Reached'; exitflag = 0;
    case 6, status = 'Memory Limit Reached'; exitflag = -4;
    case 7, status = 'Gap Limit Reached'; exitflag = 0;
    case 8, status = 'Solution Limit Reached'; exitflag = 0;
    case 9, status = 'Solution Improvement Limit Reached'; exitflag = 0;
    case 10, status = 'Restart Limit Reached'; exitflag = 0;
    case 11, status = 'Globally Optimal'; exitflag = 1;
    case 12, status = 'Infeasible'; exitflag = -1;
    case 13, status = 'Unbounded'; exitflag = -2;
    case 14, status = 'Unbounded or Infeasible'; exitflag = -2;
    otherwise, status = 'Unknown'; exitflag = -4;
end
        
        
        
        
        
        
        
