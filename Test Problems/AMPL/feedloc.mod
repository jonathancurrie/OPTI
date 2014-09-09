# feedloc.mod  LOR2-AN-90-259
# AMPL coding: S. Leyffer, January 1999.
#
# Source: Feed tray location & determination of optimum number of 
#         trays in a distillation column. Report April 1994.
#
# Number of variables:   90 (37 binary variables)  
# Number of constraints: 259 
# Objective linear
# Nonlinear constraints

param N := 12;                # max. number of trays
param M := 2;                 # number of components

set I := 1..N;
set J := 1..M;
set TOP := {N};

param F := 100;               # Feed rate
param alpha {J};              # relative volatility of component j
param xf{J};                  # feed rate of component j
param spec := 0.001;          # purity specification
param bigM := 1000;           # a large constant
param p2 := 80;               # bottom product rate

var n >= 3, <= N, integer;    # number of trays
var s{I}          binary;     # 1, if exactly i trays
var w{I}          binary;     # 1, if tray i is feed tray
var z{I}          binary;     # 1, if tray exists

var x{I,J} >= 0, <= 1;        # liquid molefraction tray i, component j
var y{I,J} >= 0, <= 1;        # vapour molefraction tray i, component j

var l  >= 0;                  # liquid mole flow rate
var v  >= 0;                  # vapour mole flow rate
var p1 >= 0, <= F;            # top product rate
var r  >= 0,  <= 5;           # reflux ratio

minimize reflux: r;

subject to

#      (1) logical constraints

       feed:  sum{i in I} w[i] = 1;

       trays: sum{i in I} s[i] = 1;
 
       defn1: sum{i in I} z[i] = n;

       defn2: sum{i in I} (i*s[i]) = n;

       order{i in I diff TOP}: z[i] - z[i+1] >= 0;

       feede: sum{i in I} (i*w[i]) <= sum{i in I} (i*s[i]);

       last{i in I}: s[i] - z[i] <= 0;

       feedi{i in I}: w[i] - z[i] <= 0;

       zto0{i in I}: z[i] <= sum{k in {i..N}} s[k];

#      (2) linear constraints

       phev1{i in I}: 0 <= sum{j in J}( y[i,j] ) + z[i] <= 2;
       phev2{i in I}: 0 <= sum{j in J}( y[i,j] ) - z[i] <= 2;
       phel1{i in I}: 0 <= sum{j in J}( x[i,j] ) + z[i] <= 2;
       phel2{i in I}: 0 <= sum{j in J}( x[i,j] ) - z[i] <= 2;

#      (3) nonlinear constraints

       eque1{i in I, j in J}: 
            ( sum{k in J}( alpha[k]*x[i,k] ) )*y[i,j] - alpha[j]*x[i,j]
            - bigM*(1 - z[i]) <= 0;

       eque2{i in I, j in J}: ( sum{k in J}( alpha[k]*x[i,k] ) )*y[i,j]
                              - alpha[j]*x[i,j] + bigM*(1 - z[i]) >= 0;

       defl: l = r*p1;

       cmbb{j in J}: p2*x[1,j] + v*y[1,j] - (l+F)*x[2,j] = 0;

       cmbt1{i in {2..N}, j in J}: l*x[i,j] + p1*y[i,j] - v*y[i-1,j] 
                                       - bigM*(1 - s[i]) <= 0;
       cmbt2{i in {2..N}, j in J}: l*x[i,j] + p1*y[i,j] - v*y[i-1,j] 
                                       + bigM*(1 - s[i]) >= 0;

       cmb1{i in {2..N-1}, j in J}: ( l + F*sum{k in {i..N}} w[k] )*x[i,j]
               + v*y[i,j] - ( l + F*sum{k in {i+1..N}} w[k] )*x[i+1,j]
               - v*y[i-1,j] - F*xf[j]*w[i] - bigM*(1 - z[i] + s[i]) <= 0;

       cmb2{i in {2..N-1}, j in J}: ( l + F*sum{k in {i..N}} w[k] )*x[i,j]
               + v*y[i,j] - ( l + F*sum{k in {i+1..N}} w[k] )*x[i+1,j]
               - v*y[i-1,j] - F*xf[j]*w[i] + bigM*(1 - z[i] + s[i]) >= 0;

       recovery{i in I}: p1*y[i,1] - bigM*(1 - s[i]) <= 0.03*xf[1];

data;

param:     alpha,    xf :=
      1    1.0       0.8
      2    5.13435   0.2 ;

