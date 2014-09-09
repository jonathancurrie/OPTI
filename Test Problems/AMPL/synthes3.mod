# synthes3.mod   OOR2-AN-11-14
# AMPL coding: S. Leyffer, December 1998.
#
# Source: Test problem 3 (Synthesis of processing system) in 
# M. Duran & I.E. Grossmann,
# "An outer approximation algorithm for a class of mixed integer nonlinear
#  programs", Mathematical Programming 36, pp. 307-339, 1986.
#
# Number of variables:   17 (8 binary variables)  
# Number of constraints: 23
# Objective nonlinear
# Nonlinear constraints

set I := 1..9;
set J := 1..8;

param u {I} default 2;

var x {i in I} >= 0, <= u[i];
var y {J}      binary;

minimize Obj:
     5*y[1] + 8*y[2] + 6*y[3] + 10*y[4] 
   + 6*y[5] + 7*y[6] + 4*y[7] + 5*y[8]
   - 10*x[1] - 15*x[2] + 15*x[3] + 80*x[4] + 25*x[5] 
   + 35*x[6] - 40*x[7] + 15*x[8] - 35*x[9]
   + exp(x[1]) + exp(0.833333*x[2]) - 65*log(x[3]+x[4]+1) 
   - 90*log(x[5]+1) - 80*log(x[6]+1) + 120;

s.t. c1:  - 1.5*log(x[5]+1) - log(x[6]+1) - x[8] <= 0;
     c2:  - log(x[3]+x[4]+1) <= 0;
     c3:  - x[1] - x[2] + x[3] + 2*x[4] + 0.8*x[5] 
          + 0.8*x[6] - 0.5*x[7] - x[8] - 2*x[9] <= 0;
     c4:  - x[1] - x[2] + 2*x[4] + 0.8*x[5] + 0.8*x[6] 
          - 2*x[7] - x[8] - 2*x[9] <= 0;
     c5:  - 2*x[4] - 0.8*x[5] - 0.8*x[6] + 2*x[7] 
          + x[8] + 2*x[9] <= 0;
     c6:  - 0.8*x[5] - 0.8*x[6] + x[8] <= 0;
     c7:  - x[4] + x[7] + x[9] <= 0;
     c8:  - 0.4*x[5] - 0.4*x[6] + 1.5*x[8] <= 0;
     c9:  0.16*x[5] + 0.16*x[6] - 1.2*x[8] <= 0;
     c10: x[3] - 0.8*x[4] <= 0;
     c11: - x[3] + 0.4*x[4] <= 0;
     c12: exp(x[1]) - 10*y[1] <= 1;
     c13: exp(0.833333*x[2]) - 10*y[2] <= 1;
     c14: x[7] - 10*y[3] <= 0;
     c15: 0.8*x[5] + 0.8*x[6] - 10*y[4] <= 0;
     c16: 2*x[4] - 2*x[7] - 2*x[9] - 10*y[5] <= 0;
     c17: x[5] - 10*y[6] <= 0;
     c18: x[6] - 10*y[7] <= 0;
     c19: x[3] + x[4] - 10*y[8] <= 0;
     c20: y[1] + y[2] = 1;
     c21: y[4] + y[5] <= 1;
     c22: - y[4] + y[6] + y[7] = 0;
     c23: y[3] - y[8] <= 0;

data;
param: u :=          # upper bounds on x[i]
            3 1
            8 1
            9 3;
