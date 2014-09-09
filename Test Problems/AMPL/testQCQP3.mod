# Test AMPL QCQP

# Number of variables: 2
# Number of constraints: 3
# Objective quadratic
# 1 quad con, 2 linear constraints

var x {1..2};

minimize Obj:
     0.5*(x[1]*x[1] + x[2]*x[2]) - 2*x[1] - 2*x[2];

s.t. c1: -x[1] + x[2] <= 2;
s.t. c2: x[1] + 3*x[2] <= 5;
s.t. c3: -x[1]*x[1] + x[2]*x[2] - 2*x[2] <= -0.5;
s.t. c4: 0 <= x[1] <= 40;
s.t. c5: 0 <= x[2];

