%% Testing Derivative Checker


%% P1
clc
fun = @(x) 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
grad = @(x)[[2*x(1)-400*x(1)*(-x(1)^2+x(2))-2];[-200*x(1)^2+200*x(2)]];
x0 = [-2; 1];        

optiDerCheck(fun,grad,x0,'Gradient')        
        

%% P1 E1 Wrong
clc
fun = @(x) 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
grad = @(x)[[2*x(1)-400*x(1)*(-x(1)^2+x(2))-2];[-20*x(1)^2+200*x(2)]]';
x0 = [-2; 1];        

optiDerCheck(fun,grad,x0,'Gradient')   


%% P2
clc
fun = @(x) (x(1)-x(2))^2 + (x(3)-1)^2 + (x(4)-1)^4 + (x(5)-1)^6;
grad = @(x)[[2*x(1)-2*x(2)];[2*x(2)-2*x(1)];[2*x(3)-2];[4*(x(4)-1)^3];[6*(x(5)-1)^5]];
nlcon = @(x) [x(1) + x(2) + x(3) + 4*x(4) - 7;
                      x(3) + 5*x(5) - 6];
nljac = @(x)[[1,1,1,4,0];[0,0,1,0,5]];
x0 = [10;7;2;-3;0.8];

optiDerCheck(fun,grad,x0,'Gradient')  
optiDerCheck(nlcon,nljac,x0,'Jacobian')  

%% P2 GE3 wrong
clc
fun = @(x) (x(1)-x(2))^2 + (x(3)-1)^2 + (x(4)-1)^4 + (x(5)-1)^6;
grad = @(x)[[2*x(1)-2*x(2)];[2*x(2)-2*x(1)];[2*x(3)-1];[4*(x(4)-1)^3];[6*(x(5)-1)^5]];
nlcon = @(x) [x(1) + x(2) + x(3) + 4*x(4) - 7;
                      x(3) + 5*x(5) - 6];
nljac = @(x)[[1,1,1,4,0];[0,0,1,0,5]];
x0 = [10;7;2;-3;0.8];

optiDerCheck(fun,grad,x0,'Gradient')  
optiDerCheck(nlcon,nljac,x0,'Jacobian')  

%% P2 J wrong
clc
fun = @(x) (x(1)-x(2))^2 + (x(3)-1)^2 + (x(4)-1)^4 + (x(5)-1)^6;
grad = @(x)[[2*x(1)-2*x(2)];[2*x(2)-2*x(1)];[2*x(3)-2];[4*(x(4)-1)^3];[6*(x(5)-1)^5]];
nlcon = @(x) [x(1) + x(2) + x(3) + 4*x(4) - 7;
                      x(3) + 5*x(5) - 6];
nljac = @(x)[[1,1.1,1,4,0];[0.1,0,1,0,6]];
x0 = [10;7;2;-3;0.8];

optiDerCheck(fun,grad,x0,'Gradient')  
optiDerCheck(nlcon,nljac,x0,'Jacobian')  

%%