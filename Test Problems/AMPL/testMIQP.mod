# Test AMPL MIQP

# Number of variables: 2
# Number of constraints: 2
# Objective quadratic
# linear constraints
# one integer constraint


var x >= 2.5 integer;
var y >= 0;
var z >= 0;

minimize Obj:
     0.5*(x*x + y*y + z*z + 2*y*z - 3*z*y) - 2*x - 3*y - z;

s.t. c1: x + y + z <= 7;
s.t. c2: 3*x - 2*y - 3*z <= 5;
s.t. c3: x - 3*y + 2*z <= 1;

