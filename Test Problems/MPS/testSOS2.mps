NAME
ROWS
 N  obj     
 L  c1      
 L  c2      
 E  c3      
COLUMNS
    x1        obj                 -1   c1                  -1
    x1        c2                   1
    x2        obj                 -2   c1                   1
    x2        c2                  -3   c3                   1
    x3        obj                 -3   c1                   1
    x3        c2                   1
    x4        obj                 -1   c1                  10
    x4        c3                 -3.5
RHS
    rhs       c1                  20   c2                  30
BOUNDS
 UP BOUND     x1                  40
 UI BOUND     x4                   3 
 LI BOUND     x4                   2
SOS
 S1 set1
    x1                  10000
    x2                  20000
    x4                  40000
 S2 set2
    x3                  30000
    x4                  40000
ENDATA