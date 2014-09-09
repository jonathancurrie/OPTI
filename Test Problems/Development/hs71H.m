function h = hs71H(x,lambda)

h = [ 2*x(4)             x(4)   x(4)   2*x(1)+x(2)+x(3);
      x(4)               0      0   x(1);
      x(4)               0      0   x(1);
      2*x(1)+x(2)+x(3)  x(1)  x(1)  0 ];
                          
h = h + lambda.ineqnonlin*-[    0      x(3)*x(4) x(2)*x(4) x(2)*x(3);
                            x(3)*x(4)     0     x(1)*x(4) x(1)*x(3);
                            x(2)*x(4) x(1)*x(4)     0     x(1)*x(2);
                            x(2)*x(3) x(1)*x(3) x(1)*x(2)     0  ];
                     
h = h + lambda.eqnonlin * diag([2 2 2 2]);