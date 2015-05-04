function ret = checkSolver(instr,err)
%DO NOT USE THIS FUNCTION - IT WILL BE REMOVED IN A FUTURE RELEASE
fprintf(2,'\n\nPLEASE UPDATE YOUR CODE TO USE optiSolver! This function will be removed in a future release.\n\n');
%Default Input Args
if(nargin < 2), err = 1; end %default to create an error
if(nargin < 1), instr = []; end
ret = [];
ret = optiSolver(instr,err);