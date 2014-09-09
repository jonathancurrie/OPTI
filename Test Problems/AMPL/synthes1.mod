# synthes1.mod   OOR2-AN-6-6
# AMPL coding: S. Leyffer, December 1998.
#
# Source: Test problem 1 (Synthesis of processing system) in 
# M. Duran & I.E. Grossmann,
# "An outer approximation algorithm for a class of mixed integer nonlinear
#  programs", Mathematical Programming 36, pp. 307-339, 1986.
#
# Number of variables:   6 (3 binary variables)  
# Number of constraints: 6
# Objective nonlinear
# Nonlinear constraints

set I := 1..3;

param u {I} default 2;

var x {i in I} >= 0, <= u[i];
var y {I}      binary;

minimize Obj:
     5*y[1] + 6*y[2] + 8*y[3] + 10*x[1] - 7*x[3] - 18*log(x[2] + 1) 
     - 19.2*log(x[1] - x[2] + 1) + 10;

s.t. c1: 0.8*log(x[2] + 1) + 0.96*log(x[1] - x[2] + 1) - 0.8*x[3] >= 0;
     c2: log(x[2] + 1) + 1.2*log(x[1] - x[2] + 1) - x[3] - 2*y[3] >= -2;
     c3: x[2] - x[1] <= 0;
     c4: x[2] - 2*y[1] <= 0;
     c5: x[1] - x[2] - 2*y[2] <= 0;
     c6: y[1] + y[2] <= 1;

data;
param: u := 
            3 1;
