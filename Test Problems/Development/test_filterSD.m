%% Test functions for filterSD [DENSE]
% J. Currie Jan 2013

% Note we have two functions:
% filtersd      - Solve a NLP with dense Jacobian using filterSD
% filtersdsp    - Solve a NLP with sparse Jacobian using filterSD

clc
clear all
%Calling form
help filtersd

%% Rosenbrock
clc
%Objective
obj = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
%Gradient
grad = @(x)[2*x(1)-400*x(1)*(x(2)-x(1)^2)-2, 200*x(2)-200*x(1)^2];

%Bounds
lb = [0;0];
ub = [5;5];

x0 = [1;0.1];
opts = []; 
opts.display = 2; 
opts.rgtol = 1e-7;

[x,f,e,i,l] = filtersd(obj,grad,x0,lb,ub,[],[],[],[],opts)

%% NLP Hock & Schittkowski #71 [1;4.743;3.8211;1.3794]
clc
%Objective & Gradient
obj = @(x) x(1)*x(4)*sum(x(1:3)) + x(3);
grad = @(x) [ x(1)*x(4) + x(4)*sum(x(1:3));
              x(1)*x(4);
              x(1)*x(4) + 1;
              x(1)*sum(x(1:3)) ];          
%Bounds
lb = ones(4,1);
ub = 5*ones(4,1);
%Nonlinear Constraints [note not transposed, using OPTI formulation]
nlcon = @(x) [ prod(x);
               sum(x.^2)];
nljac = @(x) [ prod(x)./x';
                2*x' ];    
cl = [25;40];
cu = [Inf;40];
x0 = [1 5 5 1]';
opts = []; opts.display = 2; opts.rgtol = 1e-7; 

[x,f,e,i,l] = filtersd(obj,grad,x0,lb,ub,nlcon,nljac,cl,cu,opts)

%% NLP Hock & Schittkowski #38 [1;1;1;1]
clc
%Objective & Gradient
obj = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2 + 90*(x(4)-x(3)^2)^2 + (1-x(3))^2 + ...
      10.1*(x(2)-1)^2 + 10.1*(x(4)-1)^2 + 19.8*(x(2)-1)*(x(4)-1);  
grad = @(x) [ -400*x(1)*(x(2)-x(1)^2) - 2*(1-x(1));
              200*(x(2)-x(1)^2) + 20.2*(x(2)-1) + 19.8*(x(4)-1);
              -360*x(3)*(x(4)-x(3)^2) - 2*(1-x(3));
              180*(x(4)-x(3)^2) + 20.2*(x(4)-1) + 19.8*(x(2)-1)];                          
%Bounds
lb = [-10 -10 -10 -10]';
ub = [10 10 10 10]';
x0 = [-3  -1  -3  -1]';
opts = []; opts.display = 2; opts.rgtol = 1e-7; 

[x,f,e,i,l] = filtersd(obj,grad,x0,lb,ub,[],[],[],[],opts)

%% NLP Hock & Schittkowski #51 [1;1;1;1;1]
clc
%Objective & Gradient
obj = @(x) (x(1) - x(2))^2 + (x(2) + x(3) - 2)^2 + (x(4) - 1)^2 + (x(5) - 1)^2;
grad = @(x) 2*[ x(1) - x(2);
                x(2) + x(3) - 2 - x(1) + x(2);
                x(2) + x(3) - 2;
                x(4) - 1;
                x(5) - 1 ];
%Nonlinear Constraints & Jacobian 
nlcon = @(x) [ x(1) + 3*x(2);
               x(3) + x(4) - 2*x(5);
               x(2) - x(5) ];
nljac = @(x) ([ 1  3  0  0  0;
	            0  0  1  1 -2;
	            0  1  0  0 -1 ]);
cl= [4 0 0]';
cu = [4 0 0]';
x0 = [ 2.5 0.5 2 -1 0.5 ];
opts = []; opts.display = 2; opts.rgtol = 1e-7;

[x,f,e,i,l] = filtersd(obj,grad,x0,[],[],nlcon,nljac,cl,cu,opts)

%% NLP Hock & Schittkowski #106 [f = 7.0492e+03]
clc
%Objective and Gradient
obj = @(x) x(1) + x(2) + x(3);
grad = @(x) [1 1 1 0 0 0 0 0];
%Nonlinear Constraints & Jacobian
nlcon = @(x) [1 - 0.0025*(x(4) + x(6));
              1 - 0.0025*(x(5) + x(7) - x(4));
              1 - 0.01*(x(8) - x(5));
              x(1)*x(6) - 833.33252*x(4) - 100*x(1) + 83333.333;
              x(2)*x(7) - 1250*x(5) - x(2)*x(4) + 1250*x(4);
              x(3)*x(8) - 1250000 - x(3)*x(5) + 2500*x(5)]; 
nljac = @(x) [[0,0,0,-0.0025,0,-0.0025,0,0];
              [0,0,0,0.0025,-0.0025,0,-0.0025,0];
              [0,0,0,0,0.01,0,0,-0.01];
              [x(6)-100.0,0,0,-833.33252,0,x(1),0,0];
              [0,x(7)-x(4),0,1250-x(2),-1250,0,x(2),0];
              [0,0,x(8)-x(5),0,2500-x(3),0,0,x(3)]];         
cl = [0;0;0;0;0;0];
cu = Inf(6,1);
lb = [100;1000;1000;10;10;10;10;10];
ub = [10000;10000;10000;1000;1000;1000;1000;1000];
x0 = [5000,5000,5000,200,350,150,225,425]';
opts = []; opts.display = 2; opts.rgtol = 1e-7; opts.maxtime = 0.25;

[x,f,e,i] = filtersd(obj,grad,x0,lb,ub,nlcon,nljac,cl,cu,opts)   

%% Test functions for filterSD [SPARSE]
%% NLP Hock & Schittkowski #106 [f = 7.0492e+03]
clc
%Objective and Gradient
obj = @(x) x(1) + x(2) + x(3);
grad = @(x) [1 1 1 0 0 0 0 0];
%Nonlinear Constraints & Jacobian + Structure
nlcon = @(x) [1 - 0.0025*(x(4) + x(6));
              1 - 0.0025*(x(5) + x(7) - x(4));
              1 - 0.01*(x(8) - x(5));
              x(1)*x(6) - 833.33252*x(4) - 100*x(1) + 83333.333;
              x(2)*x(7) - 1250*x(5) - x(2)*x(4) + 1250*x(4);
              x(3)*x(8) - 1250000 - x(3)*x(5) + 2500*x(5)]; 
nljac = @(x) sparse([[0,0,0,-0.0025,0,-0.0025,0,0];
                    [0,0,0,0.0025,-0.0025,0,-0.0025,0];
                    [0,0,0,0,0.01,0,0,-0.01];
                    [x(6)-100.0,0,0,-833.33252,0,x(1),0,0];
                    [0,x(7)-x(4),0,1250-x(2),-1250,0,x(2),0];
                    [0,0,x(8)-x(5),0,2500-x(3),0,0,x(3)]]);         
nljacstr =   sparse([0 0 0 1 0 1 0 0; 
                     0 0 0 1 1 0 1 0;
                     0 0 0 0 1 0 0 1;
                     1 0 0 1 0 1 0 0;
                     0 1 0 1 1 0 1 0;
                     0 0 1 0 1 0 0 1]);              
cl = [0;0;0;0;0;0];
cu = Inf(6,1);
lb = [100;1000;1000;10;10;10;10;10];
ub = [10000;10000;10000;1000;1000;1000;1000;1000];
x0 = [5000,5000,5000,200,350,150,225,425]';
opts = []; opts.display = 2; opts.rgtol = 1e-7;

%Call FilterSD [Sparse Version]
% for i = 1:50
    [x,f,e,i] = filtersdsp(obj,grad,x0,lb,ub,nlcon,nljac,nljacstr,cl,cu,opts)  
% end

%% NLP Hock & Schittkowski #51 [1;1;1;1;1]
clc
%Objective & Gradient
obj = @(x) (x(1) - x(2))^2 + (x(2) + x(3) - 2)^2 + (x(4) - 1)^2 + (x(5) - 1)^2;
grad = @(x) 2*[ x(1) - x(2);
                x(2) + x(3) - 2 - x(1) + x(2);
                x(2) + x(3) - 2;
                x(4) - 1;
                x(5) - 1 ];
%Nonlinear Constraints & Jacobian + Structure
nlcon = @(x) [ x(1) + 3*x(2);
               x(3) + x(4) - 2*x(5);
               x(2) - x(5) ];
nljac = @(x) sparse([ 1  3  0  0  0;
                      0  0  1  1 -2;
                      0  1  0  0 -1 ]);
nljacstr =   sparse([ 1  1  0  0  0;
                      0  0  1  1  1;
                      0  1  0  0  1]);
cl= [4 0 0]';
cu = [4 0 0]';
x0 = [ 2.5 0.5 2 -1 0.5 ];
opts = []; opts.display = 2; opts.rgtol = 1e-7;

%Call FilterSD [Sparse Version]
[x,f,e,i,l] = filtersdsp(obj,grad,x0,[],[],nlcon,nljac,nljacstr,cl,cu,opts)

%% NLP Hock & Schittkowski #38 [1;1;1;1]
clc
%Objective & Gradient
obj = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2 + 90*(x(4)-x(3)^2)^2 + (1-x(3))^2 + ...
      10.1*(x(2)-1)^2 + 10.1*(x(4)-1)^2 + 19.8*(x(2)-1)*(x(4)-1);  
grad = @(x) [ -400*x(1)*(x(2)-x(1)^2) - 2*(1-x(1));
              200*(x(2)-x(1)^2) + 20.2*(x(2)-1) + 19.8*(x(4)-1);
              -360*x(3)*(x(4)-x(3)^2) - 2*(1-x(3));
              180*(x(4)-x(3)^2) + 20.2*(x(4)-1) + 19.8*(x(2)-1)];                          
%Bounds
lb = [-10 -10 -10 -10]';
ub = [10 10 10 10]';
x0 = [-3  -1  -3  -1]';
opts = []; opts.display = 2; opts.rgtol = 1e-7;

[x,f,e,i,l] = filtersdsp(obj,grad,x0,lb,ub,[],[],[],[],[],opts)

%% Rosenbrock
clc
%Objective
obj = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
%Gradient
grad = @(x)[2*x(1)-400*x(1)*(x(2)-x(1)^2)-2, 200*x(2)-200*x(1)^2];

%Bounds
lb = [0;0];
ub = [5;5];

x0 = [1;0.1];
opts = []; opts.display = 2; opts.rgtol = 1e-7;

[x,f,e,stats,lambda] = filtersdsp(obj,grad,x0,lb,ub,[],[],[],[],[],opts)

%% LP [Fully Dense but run as Sparse, f = -31.4]
clc
fun = @(x) -6*x(1) - 5*x(2);
grad = @(x) [-6 -5];
nlcon = @(x) [x(1) + 4*x(2);
              6*x(1) + 4*x(2);
              2*x(1) - 5*x(2)];
nljac = @(x) sparse([1 4; 6 4; 2 -5]);
nljacstr = sparse([1 1; 1 1; 1 1]);
cl = -Inf(3,1);
cu = [16;28;6];
lb = [0;0];
ub = [10;10];

x0 = [0.1;0.5];
opts = []; opts.display = 2;
for i = 1:100
    [x,f,e,i] = filtersdsp(fun,grad,x0,lb,ub,nlcon,nljac,nljacstr,cl,cu,opts)
end

%% REMAINDER OF TESTS REQUIRE OPTI TOOLBOX

%% Large Dense Testing [use Ctrl+C to quit] [Requires OPTI Toolbox]
%Objective
% fun = @(x) sum(100*(x(1:end-1).^2-x(2:end)).^2+(x(1:end-1)-1).^2);
% grad = @(x) mklJac(fun,x);
% x0 = [zeros(1,5000) 1];
% lb = zeros(size(x0));
% ub = 5.*ones(size(x0));
% opts = []; opts.display = 2; 
% [x,fval,ef,iter] = filtersd(fun,grad,x0,lb,ub,[],[],[],[],opts);

%% Rosenbrock iterF
clc
%Objective
obj = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
%Gradient
grad = @(x) mklJac(obj,x);

%Bounds
lb = [0;0];
ub = [5;5];

x0 = [1;0.1];
opts = [];
opts.display = 2;
opts.iterfun = @optiplotlogfval;

[x,f,e,i,l] = filtersd(obj,grad,x0,lb,ub,[],[],[],[],opts)

%% AMPL IEEE System
clc
p = 'opf_014bus.nl';
prob = amplRead(p);
opts = optiset('solver','filterSD','display','final');
prob.x0 = zeros(size(prob.x0)); %make harder to solve
Opt = opti(prob,opts)
solve(Opt)

%% YALMIP Dense Failure
% clc
% n = 5;
% P = magic(n);
% Z = sdpvar(n,n,'toeplitz');
% t = sdpvar(n,n,'full');
% F = [P-Z<=t, P-Z>=-t];
% 
% opts = sdpsettings('solver','filterSD','verbose',1);
% sol = solvesdp(F,sum(sum(t)),opts);
% 
% double(Z)
% double(sum(sum(t)))

%% Johan NLP Test
% clc
% n = 10000;
% % Objective & Gradient
% fun = @(x) -sum(x);
% grad = @(x) -ones(1,n);
% 
% % Nonlinear Constraint, Jacobian & Structure 
% nlcon = @(x) x(:).^4; 
% nljac = @(x) sparse(diag(4*x(:).^3)); 
% jacstr = speye(n); 
% cl = -Inf(n,1);
% cu = ones(n,1);
% lb = -10*ones(n,1); 
% ub = 40*ones(n,1);
% 
% % Starting Guess
% x0 = zeros(n,1);
% opts = [];
% opts.display = 2;
% 
% [~,fval,exitflag] = filtersdsp(fun, grad, x0, lb, ub, nlcon, nljac, jacstr, cl, cu,opts)

