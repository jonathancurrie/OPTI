# synthes2.mod   OOR2-AN-11-14
# AMPL coding: S. Leyffer, December 1998.
#
# Source: Test problem 2(Synthesis of processing system) in 
# M. Duran & I.E. Grossmann,
# "An outer approximation algorithm for a class of mixed integer nonlinear
#  programs", Mathematical Programming 36, pp. 307-339, 1986.
#
# Number of variables:   11 (5 binary variables)  
# Number of constraints: 14
# Objective nonlinear
# Nonlinear constraints

set I := 1..6;
set J := 1..5;

param u {I} default Infinity;

var x {i in I} >= 0, <= u[i];
var y {J}      binary;

minimize Obj:
     5*y[1] + 8*y[2] + 6*y[3] + 10*y[4] + 6*y[5] 
   - 10*x[1] - 15*x[2] - 15*x[3] + 15*x[4] + 5*x[5] - 20*x[6]
   + exp(x[1]) + exp(0.833333*x[2]) - 60*log(x[4]+x[5]+1) + 140;

s.t. c1: - log(x[4]+x[5]+1) <= 0;
     c2: exp(x[1]) - 10*y[1] <= 1;
     c3: exp(0.833333*x[2]) - 10*y[2] <= 1;
     c4: 1.25*x[3] - 10*y[3] <= 0;
     c5: x[4] + x[5] - 10*y[4] <= 0;
     c6: -2*x[3] + 2*x[6] - 10*y[5] <= 0;
     c7: -x[1] - x[2] - 2*x[3] + x[4] + 2*x[6] <= 0;
     c8: -x[1] - x[2] - 0.75*x[3] + x[4] + 2*x[6] <= 0;
     c9: x[3] - x[6] <= 0;
     c10: 2*x[3] - x[4] - 2*x[6] <= 0;
     c11: -0.5*x[4] + x[5] <= 0;
     c12: -0.2*x[4] - x[5] <= 0;
     c13: y[1] + y[2] = 1;
     c14: y[4] + y[5] <= 1;

data;
param: u :=          # upper bounds on x[i]
            1 2
            2 2
            3 2
            6 3;
