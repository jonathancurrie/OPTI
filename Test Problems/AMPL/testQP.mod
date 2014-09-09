# Test AMPL QP

# Number of variables: 3
# Number of constraints: 3
# Objective quadratic
# linear constraints


param N > 0 integer, := 3;
set I := 1 .. N;

var x {i in I};

minimize Obj:
     0.5*(x[1]*x[1] + x[2]*x[2] + x[3]*x[3]) - 2*x[1] - 3*x[2] - x[3];

s.t. c1: x[1] + x[2] + x[3] <= 1;
s.t. c2: 3*x[1] - 2*x[2] - 3*x[3] <= 1;
s.t. c3: x[1] - 3*x[2] + 2*x[3] <= 1;

