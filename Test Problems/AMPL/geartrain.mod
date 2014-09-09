# geartrain1.mod   OBR2-AN-4-0
# AMPL coding: S. Leyffer, January 1999.
#
# Source: Gear Train design problem
# E. Sangren, 
# "Nonlinear Integer and Discrete Programming in Mechanical Design
#  Optimization", Trans. ASME, J. Mech. Design 112, 223-229, 1990.
#
# Number of variables:   4 (4 integer variables)  
# Number of constraints: 0
# Objective nonlinear
# Bound constraints

set I := 1..4;

var x {i in I} >= 12, <= 60, integer;

minimize GearRatio:
         ( 1.0/6.931 - (x[3]*x[2]) / (x[1]*x[4]) )^2;

