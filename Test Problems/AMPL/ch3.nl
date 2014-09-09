g3 1 1 0	# problem ch3
 3 3 1 0 3	# vars, constraints, objectives, ranges, eqns
 2 1	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 3 3 3	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 9 3	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 9 0 0 0 3	# common exprs: b,c,o,c1,o1
V3 1 0
0 2
n-1
V4 1 0
1 2
n-1
V5 1 0
2 2
n-1
V6 0 0
o0
o2
o2
n2
o0
n-1
o2
n2
v0
v3
n-1
V7 0 0
o0
o2
o2
n2
o0
n-1
o2
n2
v1
v4
n-1
V8 0 0
o0
o2
o2
n2
o0
n-1
o2
n2
v2
v5
n-1
V9 0 0
o2
o2
n2
o0
n-1
o2
n2
v0
v6
V10 0 0
o2
o2
n2
o0
n-1
o2
n2
v1
v7
V11 0 0
o2
o2
n2
o0
n-1
o2
n2
v2
v8
C0
o54
3
o2
n0.3333333333333333
v8
o2
n0.3333333333333333
v7
o2
n0.3333333333333333
v6
C1
o54
3
o2
n0.3333333333333333
v11
o2
n0.3333333333333333
v10
o2
n0.3333333333333333
v9
C2
n0
V12 1 4
0 -2
o0
v9
n1
V13 1 4
1 -2
o0
v10
n1
V14 1 4
2 -2
o0
v11
n1
O0 0
o54
3
o2
n0.5
o5
o2
n0.3333333333333333
o54
3
v3
v4
v5
n2
o2
n0.5
o5
o0
n0.3333333333333333
o2
n0.3333333333333333
o54
3
v6
v7
v8
n2
o2
n0.5
o5
o2
n0.3333333333333333
o54
3
v12
v13
v14
n2
x3
0 0.25
1 0.5
2 0.75
r
4 -0.3333333333333333
4 -1
4 1
b
3
3
3
k2
3
6
J0 3
0 0
1 0
2 0
J1 3
0 -0.6666666666666666
1 -0.6666666666666666
2 -0.6666666666666666
J2 3
0 0.6666666666666666
1 0.6666666666666666
2 0.6666666666666666
G0 3
0 0
1 0
2 0
