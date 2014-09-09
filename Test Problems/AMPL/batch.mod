# batch.mod   OOR2-AN-46-73
# AMPL coding: S. Leyffer, January 1999.
#
# Source: Source: Optimal Design of Multiproduct Batch Plant
#         G.R. Kocis & I.E. Grossmann,  "Global Optimization 
#         of Nonconvex Mixed Integer Nonlinear Programmming
#         (MINLP) problems in Process Synthesis", Indust. Engng. Chem. Res.,
#         No. 27, pp 1407--1421, 1988.
#
# Number of variables:   46 (24 binary variables)  
# Number of constraints: 73 
# Objective nonlinear
# Nonlinear constraints

param m := 6;            # Number of stages
param n := 5;            # Number of products
param k := 4;            # Max number of parallel usints per stage

set M := 1..m;
set N := 1..n;
set K := 1..k;
set KxM := (K cross M);
set NxM within (N cross M);

param alpha{M} default 250;  # cost coefficient
param beta{M}  default 0.6;  # cost coefficient
param VL{M}    default 300;  # lower bound on volume
param VU{M}    default 3000; # upper bound on volume
param S{NxM};                # size factor (L/kg)
param t{NxM};                # processing time (h)
param Q{N};                  # production rate (kg)
param tLO{N};                # lower bounds on cycle time
param tUP{N};                # upper bounds on cycle time
param bLO{N};                # lower bounds on product rate
param bUP{N};                # upper bounds on product rate
param H := 6000;             # horizon time (h)

var npu {M} >= 0, <= log(k); # log of number of parallel units at stage j
var NPU {M} >= 0, <= k;      # number of parallel units at stage j

var v {j in M} >= log(VL[j]), <= log(VU[j]); # log of volume at stage j
var b {i in N} >= bLO[i], <= bUP[i];         # log of batch size of batch i
var tL{i in N} >= tLO[i], <= tUP[i];         # log of cycle time of batch i

var y {KxM} binary;  # models number of parallel units

minimize cost:
      sum {j in M} ( alpha[j] * exp(npu[j] + beta[j]*v[j]) );

subject to
      volm {i in N, j in M}:  v[j] >= log(S[i,j]) + b[i];
      
      ctim {i in N, j in M}: npu[j] + tL[i] >= log(t[i,j]);

      horizon: sum {i in N} ( Q[i]*exp(tL[i] - b[i]) ) <= H;

      npar {j in M}: npu[j] = sum {l in K} ( log(l) * y[l,j] );

      SOS1 {j in M}: sum {l in K} ( y[l,j] ) = 1;

data;

param:  NxM:     S,     t :=
        1 1     7.9    6.4
        2 1     0.7    6.8
        3 1     0.7    1.0
        4 1     4.7    3.2
        5 1     1.2    2.1
        1 2     2.0    4.7
        2 2     0.8    6.4
        3 2     2.6    6.3
        4 2     2.3    3.0
        5 2     3.6    2.5
        1 3     5.2    8.3
        2 3     0.9    6.5
        3 3     1.6    5.4
        4 3     1.6    3.5
        5 3     2.4    4.2
        1 4     4.9    3.9
        2 4     3.4    4.4
        3 4     3.6    11.9
        4 4     2.7    3.3
        5 4     4.5    3.6
        1 5     6.1    2.1
        2 5     2.1    2.3
        3 5     3.2    5.7
        4 5     1.2    2.8
        5 5     1.6    3.7
        1 6     4.2    1.2
        2 6     2.5    3.2
        3 6     2.9    6.2
        4 6     2.5    3.4
        5 6     2.1    2.2;

param:       Q  := 
       1  250000.0    
       2  150000.0    
       3  180000.0    
       4  160000.0    
       5  120000.0;

param:      tLO :=
       1  0.729961     
       2  0.530628 
       3  1.09024 
       4  -0.133531 
       5  0.0487901 ;

param:      tUP :=
       1  2.11626            
       2  1.91626
       3  2.47654
       4  1.25276
       5  1.43508;

param:      bLO :=
       1  4.45966    
       2  3.74950 
       3  4.49144
       4  3.14988
       5  3.04452;

param:      bUP :=
       1  397.747     
       2  882.353
       3  833.333
       4  638.298
       5  666.667;
