# MIQP3
var x >= 0, <= 10 integer;
var y >= 0, <= 10;

minimize Obj:
     0.5*(x*x + 2*y*y - 2*x*y) - 2*x - 6*y;

s.t. c1: x + y <= 3;
s.t. c2: -x + 2*y <= 3.5;


