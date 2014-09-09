# trimloss.mod   LOR2-AN-V-V
# AMPL coding: S. Leyffer, January 1999.
#
# Source: I. Harjunkoski, T. Westerlund, R. P\"{o}rn and H. Skrifvars
# "Different transformations for solving non--convex trim loss problems 
# by MINLP", European Journal of Operational Research 105 (1998) 594-603.
#
# Convexified version of trimlon.mod (using square root transformation).
#
# Convex MINLP arising from
# trim loss minimization in the paper industry. The problem is to 
# produce a set of product paper rolls from raw paper rolls such that 
# a cost function including the trim loss and the overall production cost
# is minimized.
#
# This formulation requires any of the following data files:
#
#                  variables   constraints   comments
#    ----------------------------------------------------------------
#    trimloss2.dat
#    trimloss4.dat
#    trimloss5.dat
#    trimloss6.dat
#    trimloss7.dat
#    trimloss12.dat
#
# Number of variables:   variable
# Number of constraints: variable
# Objective linear
# Quadratic constraints

set I;
set J;

param L{J} integer;		# max. multiple of cutting pattern
param K{I} integer;		# max. no. of products i in cut

param c{J};                     # cost of raw material
param C{J};                     # cost of change-over of knives
param b{I};                     # width of product roll-type i
param Bmax;                     # width of raw paper roll
param Delta;                    # tolerance for width
param Nmax;                     # max number of products in cut
param nord{I};                  # number of orders of product type i

var y{J} binary;    		# = 1, if cutting pattern present
var beta{j in J,1..L[j]} binary;# SOS for max. multiple of cutting pattern
var z{i in I,J,1..K[i]} binary;	# SOS for max. no. of products i in cut
var M{J} >= 1;      		# square of no. repeats of pattern j
var N{I,J} >= 1; 		# square of no. products i produced in cut j

minimize cost: sum{j in J} ( c[j]*( sum{l in 1..L[j]} (beta[j,l]*l) ) 
                             + C[j]*y[j] );
         
subject to

   max_width{j in J}:   sum{i in I} ( b[i]* (sum{k in 1..K[i]} z[i,j,k]*k) ) 
			<= Bmax;         

   min_width{j in J}:   sum{i in I} ( b[i]* (sum{k in 1..K[i]} z[i,j,k]*k) )
			>= Bmax - Delta;    

   max_n_sum{j in J}:   sum{i in I} ( sum{k in 1..K[i]} z[i,j,k]*k ) <= Nmax;

   cut_exist{j in J}:   y[j] - sum{l in 1..L[j]} (beta[j,l]*l) <= 0;

   no_cut{j in J}:      sum{l in 1..L[j]} (beta[j,l]*l) - L[j]*y[j] <= 0;

   def_M{j in J}:	M[j] = 1 + sum{l in 1..L[j]}( beta[j,l]*l*(l+2) );
  
   SOS_beta{j in J}:	sum{l in 1..L[j]} beta[j,l] <= 1;

   def_N{i in I, j in J}:
			N[i,j] = 1 + sum{k in 1..K[i]}( z[i,j,k]*k*(k+2) );

   SOS_z{i in I, j in J}:
			sum{k in 1..K[i]} z[i,j,k] <= 1;

   order{i in I}:	nord[i] + card(J) - sum{j in J} sqrt(M[j]*N[i,j])
			+ sum{j in J}( sum{l in 1..L[j]}(beta[j,l]*l)
					+ sum{k in 1..K[i]}(z[i,j,k]*k) ) <= 0;
