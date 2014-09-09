# c-sched.mod   OLR2-AN-v-v
# AMPL coding: S. Leyffer, University of Dundee, March 2000.
#
# Source: V. Jain & I.E. Grossmann, "Cyclic Scheduling of Continuous
#         Parallel Units with Decaying Performance", AIChE Journal, 
#         44, 1623-1636.
#
# Formulation with binary variables
#
# Two data files: c-sched1.dat	= Example 1 (3 feeds, single furnace)
#		  c-sched2.dat  = Example 2 (7 feeds, 4 furnaces)

# ... index sets
set feeds;				# set of feeds
set furnaces;				# set of furnaces

# ... parameters (constants)
param K;				# upper bound on No. subcycles
param U;				# upper bound on processing time
param epsi := 0.01;			# small const to avoid 0/0
param tau {feeds,furnaces};		# changeover time [days]
param D {feeds,furnaces};		# processing rate [tons/day]
param a {feeds,furnaces};		# conversion parameter [1/day]
param b {feeds,furnaces};		# conversion parameter [1/day]
param c {feeds,furnaces};		# conversion parameter [1/day]
param P {feeds,furnaces};		# price parameter [$/ton]
param Cs {feeds,furnaces};		# setup/cleaning cost [$]
param Flo {feeds};			# lower bnd on flow rate
param Fup {feeds};			# upper bnd on flow rate

# ... derived parameters
param Cc {i in feeds, l in furnaces} := P[i,l]*D[i,l]*c[i,l];
param Cp {i in feeds, l in furnaces} := P[i,l]*D[i,l]*a[i,l] / b[i,l];
param yk {k in 0..K} := if k == 0 then epsi else k;

# ... problem variables
var t {feeds,furnaces} >= 0;		# process time of feed in furnace
var n {feeds,furnaces} >= 0, <= K, 
      := 1.0;				# number of subcycles of feed in furnace
var F {i in feeds} >= Flo[i], <= Fup[i];# rate of arrival of feed i
var S {feeds} >= 0;			# extra amount of feed processed above min
var dt {feeds,furnaces} >= 0;		# time devoted to feed in furnace
var Tcycle >= 0;			# common cycle time for all furnaces
var y {feeds,furnaces,0..K} >= 0, <= 1,
		            binary;	# SOS to model n[i,l]

# ... objective is to maximize profit/cycle-time
maximize profit: sum{i in feeds, l in furnaces}
                    (   Cc[i,l]*t[i,l] - Cs[i,l]*n[i,l]
                      + Cp[i,l]*n[i,l]*(1 - exp(-b[i,l]*t[i,l]/n[i,l]) )
                    ) / Tcycle;
                      
subject to

   # ... mass balance equations (8) & (9) 
   massbal_1{i in feeds}: Flo[i]*Tcycle + S[i] = sum{l in furnaces}( D[i,l]*t[i,l] );
   massbal_2{i in feeds}: S[i] <= ( Fup[i] - Flo[i] )*Tcycle;

   # ... integrality constraints (10) & (11)
   integ_1{i in feeds, l in furnaces}: n[i,l] = sum{k in 0..K} yk[k] * y[i,l,k];
   integ_2{i in feeds, l in furnaces}: 1      = sum{k in 0..K} y[i,l,k];

   # ... timing constraint (13): relate total time of feed to processing & clean-up
   time_1{i in feeds, l in furnaces}: dt[i,l] = n[i,l]*tau[i,l] + t[i,l];

   # ... timing constraint (14): total time less than cycle time
   time_2{l in furnaces}: sum{i in feeds} dt[i,l] <= Tcycle;

   # ... timing constraint (15): t[i,l] is zero if number of subcycles is zero
   time_3{i in feeds, l in furnaces}: t[i,l] <= U*(1 - y[i,l,0]);

   # .. extra constraints
   extra{i in feeds: Flo[i] > 0}: sum{l in furnaces} n[i,l] >= 1;
   