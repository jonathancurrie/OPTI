# lbti-00.mod  
#
# Source: Mark A. Abramson "Mixed Variable Optimization of a
# Load-Bearing Thermal Insulation System Using a Filter Pattern Search
# Algorithm", Optimization and Engineering, 5, 157-177, 2004.  
#
# MINLP version of the model implements the model described in
# K. Abhishek, S. Leyffer, and J.T. Linderoth, Modeling without Categorical
# Variables: A Mixed Integer Nonlinear Program for the Optimization of Thermal
# Insulation Systems, Preprint ANL/MCS-Pnnnn-mmmm. The equations numbers (i.j)
# refer to the equation numbers in that report.
#
# This formulation requires the following data files:
#    nylon.dat
#    teflon.dat
#    steel.dat
#    carbon-steel.dat
#    aluminum.dat
#    epoxy-normal.dat
#    epoxy-plane.dat
#
# Difference to lbti-01.mod: DISCONTINUOUS OBJECTIVE COEFFICIENTS
#
# Objective nonlinear
# Some Constraints are nonlinear
# Version 0.1 (Oct. 27, 2006)
# Version 0.2 (Jan. 29, 2007) ... Sven corrected bug in Defn of v[i,j], w[i,j]

# ... sets 
set materials;
set epoxy within materials;              # ... handle epoxy material specially

# ... model parameters (constants, see Table 1)
param DELTA >= 0, <= 5.0, default 0.05;  # max % expansion allowed.
param F >= 0, default 250000;            # max allowed force on the system.
param L >= 0, default 100;               # length of the system.
param M >= 0, default 10;                # max mass of the system.
param N >= 1, integer, default 10;       # max number of intercepts.
param SIGMAX >= 0, default 100000;
param SIGMIN >= 0, default 1000;
param TCOLD >= 0, default 4.2;           # temperature at the coldest end.
param THOT >= 0, default 300;            # temperature at the hottest end.
param AMAX >= 0, default F/SIGMIN;       # maximum area of insulator
param AMIN >= 0, default F/SIGMAX;       # minimum area of insulator
param D{materials} integer, default 20;  # number of discretization points of data splines
param EPS > 0, default 1E-2;             # small constant 

# ... model data (lookup tables/splines, see Table 1)
param RHO{materials};                    # materials densities 
param T{j in materials, 1..D[j]};        # temperature data (discretization points)
param K{j in materials, 1..D[j]};        # thermal conductivity data
param SIGMA{j in materials, 1..D[j]};    # tensile yield strength data
param E{j in materials, 1..D[j]};        # thermal expansion of insulator

# ... integral of k() in each discretization interval with trapezoidal rule
param AK{j in materials, r in 1..D[j]-1} 
      := 0.5*(K[j, r] + K[j, r+1])*(T[j,r+1] - T[j,r]); 
# ... integral of e().k() in each discretization interval with trapezoidal rule
param AEK{j in materials, r in 1..D[j]-1} 
      := 0.5*(E[j, r] + E[j, r+1])*(T[j,r+1] - T[j,r]); 

# ... binary/integer VARIABLES
var n integer, >= 0, <= N;          # number of intercepts to use.
var y{1..N+1} binary;	            # indicates if intercept exists
var z{1..N+1, materials} binary;    # indicate material of insulator
				    # used between intercepts i-1 and i  

# ... continuous VARIABLES
var x{1..N+1} := L/(N+1),  >= 0, <= L;     # thickness of each insulator
var t{0..N+1} >= TCOLD, <= THOT;           # temperature of each intercept.
var a{1..N+1} := AMIN , >= 0, <= AMAX;     # area of each insulator 
var u{1..N+1} >= 0;                        # variables to do the transformation 
                                           # u=w.x/a. 
var q{1..N+1} >= 0;                             # heat flow from intercept i to i-1
var power >= 0;                            # cooling power 

# ... defined variable ... approx. integral from t[i-1] to t[i] of AK (3.30)
var v{i in 1..N+1, j in materials}
    = sum{r in 1..D[j]-1}
	if (T[j,r] >= t[i-1] and T[j,r+1] <= t[i]) then
	   AK[j, r]
	else if (T[j,r] < t[i-1] and T[j,r+1] > t[i-1] and T[j,r+1] <= t[i]) then
	   (  K[j,r] + (K[j,r+1]-K[j,r])/(T[j,r+1]-T[j,r])*(t[i-1]-T[j,r])
	    + K[j,r+1] )*( T[j,r+1]-t[i-1] )/2
        else if (T[j,r] < t[i] and T[j,r+1] > t[i] and T[j,r] >= t[i-1]) then
	   (  K[j,r] + (K[j,r+1]-K[j,r])/(T[j,r+1]-T[j,r])*(t[i]-T[j,r])
            + K[j,r] )*( t[i] - T[j,r] )/2 
        else if (T[j,r] < t[i-1]  and T[j,r+1] > t[i]) then
           (  K[j,r] + (K[j,r+1]-K[j,r])/(T[j,r+1]-T[j,r])*(t[i-1]-T[j,r])
            + K[j,r] + (K[j,r+1]-K[j,r])/(T[j,r+1]-T[j,r])*(t[i  ]-T[j,r])
           )*( t[i] - t[i-1] )/2
    ;

# ... defined variable ... approx. integral from t[i-1] to t[i] of AEK (3.31)
var w{i in 1..N+1, j in materials}
    = sum{r in 1..D[j]-1} 
	if (T[j,r] >= t[i-1] and T[j,r+1] <= t[i]) then
	   AEK[j, r]
	else if (T[j,r] < t[i-1] and T[j,r+1] > t[i-1] and T[j,r+1] <= t[i]) then
	   (  E[j,r] + (E[j,r+1]-E[j,r])/(T[j,r+1]-T[j,r])*(t[i-1]-T[j,r])
	    + E[j,r+1] )*( T[j,r+1]-t[i-1] )/2
        else if (T[j,r] < t[i] and T[j,r+1] > t[i] and T[j,r] >= t[i-1]) then
	   (  E[j,r] + (E[j,r+1]-E[j,r])/(T[j,r+1]-T[j,r])*(t[i]-T[j,r])
            + E[j,r] )*( t[i] - T[j,r] )/2 
        else if (T[j,r] < t[i-1]  and T[j,r+1] > t[i]) then
           (  E[j,r] + (E[j,r+1]-E[j,r])/(T[j,r+1]-T[j,r])*(t[i-1]-T[j,r])
            + E[j,r] + (E[j,r+1]-E[j,r])/(T[j,r+1]-T[j,r])*(t[i  ]-T[j,r])
           )*( t[i] - t[i-1] )/2
    ;

# ... defined discontinuous thermal-dynamic cycle efficiency (2.1)
var C{i in 0..N+1} = if      ( t[i] <=  4.2 ) then 5.0
                     else if ( t[i] <  71   ) then 4.0
                     else                          2.5 ;

minimize objf: power;

subject to

	# ... definition of objective (paper, see results/lbti-00.out-1)
	## DefP: power = sum{i in 1..N} y[i] * C[i] * ( THOT/t[i] - 1 ) * ( q[i+1] - q[i] );
	## DiffQ{i in 1..N}: q[i+1] - q[i] >= 0;
	# ... definition of objective (matlab, see results/lbti-00.out-2)
	DefP: power =   sum{i in 0..N} y[i+1] * C[i] * (THOT/t[i] - 1) * q[i+1] 
                      - sum{i in 1..N+1} y[i] * C[i] * (THOT/t[i] - 1) * q[i];

	# ... definition of number of intercepts (3.1)
	Defn: sum{i in 1..N+1} y[i] = n + 1;

	# ... ensure intercepts are contiguous (3.2)
	Ordy{i in 1..N}: y[i+1] <= y[i]; 

	# ... force x[i]=0 for i>n+1 (3.4)
	xUpper{i in 1..N+1}: sum{j in i..N+1} x[j] <= L*y[i];

	# ... pick exactly one material for each intercept (3.13)
	ChooseMaterial1{i in 1..N+1}: sum{j in materials} z[i,j] = y[i];

	# ... max total length of insulators (3.14)
	length: sum{i in 1..N+1} x[i] = L;

	# ... total mass constraint (3.16)
	mass: sum{i in 1..N+1, j in materials} RHO[j]*z[i,j]*a[i]*x[i] <= M;

	# ... force x[i]>L*EPS for i<=n+1 (3.18)
	xLower{i in 1..N+1}: x[i] >= EPS*L*y[i];

	# ... set temperature at boundary (3.19)
	Tcold:           t[0]   = TCOLD;
	THot{i in 1..N}: t[i]  >= THOT*(1 - y[i+1]);	
	TN:              t[N+1] = THOT;

        # ... force temperatures apart for existing layers (3.20)
	TempDiff{i in 1..N+1}: t[i] - t[i-1] >= EPS*y[i];

	# ... upper/lower bound on area of existing intercepts (3.22)
	#MaxAreaConstr{i in 1..N+1}: sum{j in i..N+1} a[j] <= AMAX*y[i];
	#MinAreaConstr{i in 1..N+1}: sum{j in i..N+1} a[j] >= AMIN*y[i];
	MaxAreaConstr1{i in 1..N+1}: a[i] <= AMAX*y[i];
	MinAreaConstr1{i in 1..N+1}: a[i] >= AMIN*y[i];

	# ... area stress constraints for all (monotone) materials (3.26)
	#     (force per unit area for the material less than allowable stress)
	areaStressAll{i in 1..N+1, j in materials}:
	a[i]*sum{r in 1..D[j]-1}(
		if (t[i] >= T[j,r] and t[i] <= T[j,r+1]) then 
		   SIGMA[j,r] + ((SIGMA[j,r+1] - SIGMA[j,r])/(T[j,r+1] - T[j,r]))
                                *(t[i] - T[j,r])) 
        >= F*z[i,j] ;

	# ... area stress constraints for expoxys (3.27)
	#     (force per unit area for the material less than allowable stress)
	areaStressEpoxyN{i in 1..N+1, j in epoxy}: 
	a[i]*sum{r in 1..D[j]-1}( 
		if (t[i-1] >= T[j,r] and t[i-1] <= T[j,r+1]) then 
		   SIGMA[j,r] + ((SIGMA[j,r+1] -SIGMA[j,r])/(T[j,r+1] - T[j,r]))
                                *(t[i-1] - T[j,r])) 
	>= F*z[i,j];

	# ... heat flow is governed by Fourier's law (3.32) 
	# changed = to <= (debug 23th feb)
	HeatFlow{i in 1..N+1}: 
	a[i] * sum{j in materials} v[i,j]*z[i,j] <= x[i]*q[i] ;  

	# ... thermal contraction $\Delta x_i / x_i$ (3.33)
	# changed = to <= (debug 23th feb)
	ThermalContraction{i in 1..N+1}: 
	sum{j in materials} w[i,j]*z[i,j] <= u[i]*sum{j in materials} v[i,j]*z[i,j];

	# ... thermal expansion/contraction constraint (3.34)
	ThermalExpansion: sum{i in 1..N+1} u[i]*x[i] <= L*DELTA/100;	      
