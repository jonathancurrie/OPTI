# c-reload.mod   OOR2-AN-342-284
# AMPL coding: S. Leyffer, University of Dundee, April 2000.
#
# Source: 
# "Optimization of a Nuclear Reactor Core Reload Pattern Using
#  Nonlinear Optimization and Search Heuristics", A.J. Quist,
#  R. van Geemert, J.E. Hoogenboom, T. Illes, E. de Klerk, T. Terlaki
#  Delft University.
#
# Data files: c-reload-1.dat

set I;				# ... set of nodes
set L ordered;			# ... set of batches
set M;				# ... set of bundles in batches
set T ordered;			# ... set of time steps
set ND := 1..2;			# ... number of diagonal pairs
set EoC := { last(T) };		# ... end of cycle definition

# ... scalar parameters
param conspw;			# ... conspw = Pc in the article
param kfresh;
param alfa;
param flim;
param cytime;			# ... Delta t = cytime/(card(T)-1) = 70 here

# ... array parameters
param d1{ND} integer;		# ... diagonal pairs
param d2{ND} integer;		# ... x[d1[i],l,m] = x[d2[i],l,m]
param v{I} default 1.0;
param G{I,I} default 0.0;	# ... matrix G[i,j]

# ... variables
var r {I,T}    >= 0, := kfresh/14;	# ... phi in paper	
var x {I,L,M} binary, := 0.5;
var kinf {I,T} >= 0, := kfresh;
var keff {T}   >= 0;
var phi  {I,T} >= 0;

### OBJECTIVE FUNCTION ###
maximize goal: sum{e in EoC} keff[e];

subject to

  # ... every bundle is assigned in 1 node (or 2 halves)
  sumi{l in L, m in M}: sum{i in I} v[i]*x[i,l,m] = 1 ;

  # ... every node contains 1 bundle
  sumlm{ i in I}: sum{l in L, m in M} x[i,l,m] = 1 ;

  # ... standard tri-linear reloading equation
  plac{i in I, e in EoC}: 
    kinf[i,1] = sum{l in L: l > 1}
                ( sum{m in M} x[i,l,m]*( sum{j in I} v[j]*x[j,l-1,m]*kinf[j,e] ) )
                + ( sum{m in M} x[i,1,m] ) * kfresh;

  # ... standard kernel equation definition
  kern{i in I, t in T}: sum{j in I} (G[i,j]*kinf[j,t]*r[j,t]) = r[i,t]*keff[t];

  # ... burnup described as decay of kinf
  kinff{i in I, t in T diff EoC}: 
       kinf[i,t+1] = kinf[i,t]
                     - alfa*conspw*cytime/(card(T)-1) * kinf[i,t]*r[i,t];

  # ... normalize the power density to 1
  cpow{t in T}: sum{i in I} v[i]*kinf[i,t]*r[i,t] = 1;

  # ... power peaking constraint
  peak{i in I, t in T}: kinf[i,t]*r[i,t] <= flim / (sum{j in I} v[j]);

  # ...fills adjacent diagonal elements with the same bundle
  diag{n in ND, l in L, m in M}: x[d1[n],l,m] = x[d2[n],l,m];
