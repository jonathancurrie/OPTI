# Test AMPL NLP with inf

# Number of variables: 1
# Number of constraints: 1
# Objective nonlinear
# linear constraints


param N > 0 integer, := 1;
set I := 1 .. N;

var x {i in I};

maximize Obj:
     1/x[1];

s.t. c1: 0 <= x[1] <= 1;

