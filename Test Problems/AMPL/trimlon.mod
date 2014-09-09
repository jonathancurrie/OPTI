# trimlon.mod   LQR2-AN-V-V
# AMPL coding: S. Leyffer, January 1999.
#
# Source: I. Harjunkoski, T. Westerlund, R. P\"{o}rn and H. Skrifvars
# "Different transformations for solving non--convex trim loss problems 
# by MINLP", European Journal of Operational Research 105 (1998) 594-603.
#
# Nonconvex MINLP arising from
# trim loss minimization in the paper industry. The problem is to 
# produce a set of product paper rolls from raw paper rolls such that 
# a cost function including the trim loss and the overall production cost
# is minimized.
#
# This formulation requires any of the following data files:
#
#                  variables   constraints   comments
#    ----------------------------------------------------------------
#    trimlon2.dat
#    trimlon4.dat
#    trimlon5.dat
#    trimlon6.dat
#    trimlon7.dat
#    trimlon12.dat
#
# Number of variables:   variable
# Number of constraints: variable
# Objective linear
# Quadratic constraints

set I;
set J;

param c{J};                     # cost of raw material
param C{J};                     # cost of change-over of knives
param b{I};                     # width of product roll-type i
param Bmax;                     # width of raw paper roll
param Delta;                    # tolerance for width
param Nmax;                     # max number of products in cut
param nord{I};                  # number of orders of product type i
param M;                        # max number of repeats of any pattern

var y{J} >= 0, <= 1, binary;             # = 1, if cutting pattern present
var m{J} >= 0, <= M, integer, := 1;      # number of repeats of pattern j
var n{I,J} >= 0, <= Nmax, integer, := 1; # number of products i produced in cut j

minimize cost: sum{j in J} ( c[j]*m[j] + C[j]*y[j] );
         
subject to
   max_width{j in J}:   sum{i in I} ( b[i]* n[i,j] ) - Bmax <= 0;         
   min_width{j in J}: - sum{i in I} ( b[i]* n[i,j] ) + Bmax - Delta <= 0;    
   max_n_sum{j in J}:   sum{i in I} ( n[i,j] ) - Nmax <= 0;
   cut_exist{j in J}:   y[j] - m[j] <= 0;
   no_cut{j in J}:      m[j] - M*y[j] <= 0;
   min_order{i in I}:   nord[i] - sum{j in J} ( m[j] * n[i,j] ) <= 0;
