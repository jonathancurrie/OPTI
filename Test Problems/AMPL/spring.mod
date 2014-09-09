# spring.mod   OOR2-AN-17-10
# AMPL coding: S. Leyffer, January 1999.
#
# Source: E. Sangren, 
# "Nonlinear Integer and Discrete Programming in Mechanical Design
#  Optimization", Trans. ASME, J. Mech. Design 112, 223-229, 1990.
#
#   Coil compression spring design problem, finds minimum volume 
#   of wire for the production of a coil compression spring.
#   Using special ordered sets.
#
# Number of variables:   17 (11 binary variables = 1 SOS1)  
# Number of constraints: 10
# Objective nonlinear
# Nonlinear constraints

set I := 1..11;

param b{I};                     # allowable wire diameters (in)
param Pload := 300;		# preload (lb)
param Pmax  := 1000;		# max working load (lb)
param delm  := 6;		# max deflection (in)
param delw  := 1.25;		# delection from preload (in)
param lmax  := 14;		# max free length of spring (in)
param Dmax  := 3;		# max coil diameter (in)
param S     := 189000;		# max shear stress
param G     := 11.5E6;		# shear modulus of material
param dmin  := 0.2;		# min wire diameter
param pi    := 3.141592654;	# pi

var y{I} binary;                # SOS1 for wire diameter
var D >= 2*dmin;                # coil diameter
var d >= dmin, <= 0.5;          # wire diameter
var N >= 1,   <= 100, integer;  # number of coils
var C >= 1.1;			# = D / d (D >= d)
var K;				# = (4C-1)/(4C-4)+0.615/C
var del >= 0;   		# = 8 * (N * D^3) / (G*d^4)

minimize material:
         pi * D * d^2 * (N+2)/4;

subject to
   Def_C:    C = D / d;
   Def_K:    K = (4*C - 1) / (4*C - 4) + 0.615 / C;
   shear:    S >= 8*Pmax * K * D / (pi * d^3);
   Def_del:  del = 8 * (N * D^3) / (G * d^4);
   free_l:   lmax >= Pmax*del + 1.05*(N + 2)*d;
   diametr:  Dmax >= D + d;
   deflect1: delm >= Pload * del;
   deflect2: delw <= (Pmax - Pload) * del;
   Def_d:    d = sum {i in I} ( b[i] * y[i] );	## Definition of SOS1 or SOS3
   SOS1:     1 = sum {i in I} (        y[i] );	## convexity constraint
         
data;

param:       b  :=      # allowable diameters
        1    0.207
        2    0.225
        3    0.244
        4    0.263
        5    0.283
        6    0.307
        7    0.331   
        8    0.362
        9    0.394
       10    0.4375
       11    0.500 ;

suffix sosno  IN;
let {i in I} y[i].sosno := 1;
suffix ref IN;
let {i in I} y[i].ref   := b[i];


