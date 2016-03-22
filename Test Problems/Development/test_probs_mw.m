%% LP1
clc
%Objective & Constraints
f = -[6 5]';
A = ([1,4; 6,4; 2, -5]); 
b = [16;28;6];    
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'bounds',[0;0],[10;10])
[x,fval,exitflag,info] = solve(Opt) 
%Plot
plot(Opt)
%Check Solution
[ok,msg] = checkSol(Opt)

%%
[x,fval,exitflag,info,lambda,Opt] = opti_linprog(f,A,b,[],[],[0;0],[10;10])

%% LP2
clc
%Objective & Constraints
f = -[1 2 3]';
A = [-1,1,1; 1,-3,1];
b = [20,30]';
%Setup Options
opts = optiset('solver','clp');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'eq',[1 1 1],40,'bounds',[0 0 0]',[40 inf inf]','options',opts)
[x,fval,exitflag,info] = solve(Opt)

%%
[x,fval,exitflag,info,lambda,Opt] = opti_linprog(f,A,b,[1 1 1],40,[0;0;0],[40;inf;inf],[],opts)

%% LP3
clc
%Objective & Constraints
f = [8,1]';
A = [-1,-2;1,-4;3,-1;1,5;-1,1;-1,0;0,-1]; 
b = [-4,2,21,39,3,0,0]';
%Setup Options
opts = optiset('solver','scip');
%Build & Solve
x0 = [9,6]';
Opt = opti('grad',f,'ineq',A,b,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)
%Plot
plot(Opt,10)

%%
[x,fval,exitflag,info,lambda,Opt] = opti_linprog(f,A,b,[],[],[],[],[],opts)

%% MILP1
clc
%Objective & Constraints
f = -[6 5]';
A = [1,4; 6,4; 2, -5]; 
b = [16;28;6];  
%Setup Options
opts = optiset('solver','glpk');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'bounds',[0;0],[10;10],'int','II','options',opts)
[x,fval,exitflag,info] = solve(Opt)
%Plot
plot(Opt)

%%
[x,fval,exitflag,info,Opt] = opti_mintprog(f,A,b,[],[],[0;0],[10;10],'II',[],opts)
[x,fval,exitflag,info,Opt] = opti_intlinprog(f,[1 2],A,b,[],[],[0;0],[10;10],opts)

%% MILP2
clc
%Objective & Constraints
f = -[1 2 3 1]'; 
A = [-1 1 1 10; 1 -3 1 0]; 
b = [20;30];  
Aeq = [0 1 0 -3.5];
beq = 0;
%Setup Options
opts = optiset('solver','cbc');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'eq',Aeq,beq,'bounds',[0 0 0 2]',[40 inf inf 3]','int','CCCI','options',opts)
[x,fval,exitflag,info] = solve(Opt)

%%
[x,fval,exitflag,info,Opt] = opti_mintprog(f,A,b,Aeq,beq,[0;0;0;2],[40;inf;inf;3],'CCCI',[],opts)
[x,fval,exitflag,info,Opt] = opti_intlinprog(f,4,A,b,Aeq,beq,[0;0;0;2],[40;inf;inf;3],opts)

%% MILP3 Infeasible
clc
%Objective & Constraints
f = -[6 5]';
A = [4,1; 1.5,1; -2, 1; -0.2, 1]; 
b = [5.3;4.5;-2.5; 1.5];  
%Setup Options
opts = optiset('solver','lp_solve');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'bounds',[0;0],[10;10],'int','II','options',opts)
[x,fval,exitflag,info] = solve(Opt)
%Plot
plot(Opt)

%%
[x,fval,exitflag,info,Opt] = opti_mintprog(f,A,b,[],[],[0;0],[10;10],'II',[],opts)
[x,fval,exitflag,info,Opt] = opti_intlinprog(f,[1 2],A,b,[],[],[0;0],[10;10],opts)

%% MILP4 SOS
clc
%Objective & Constraints
f = [-1 -1 -3 -2 -2]';
A = [-1 -1 1 1 0;
      1 0 1 -3 0];
b = [30;30];
lb = zeros(5,1);
ub = [40;1;inf;inf;1];
sos = '1';
sosind = [1:5]';
soswt = [1:5]';
%Setup Options
opts = optiset('solver','cbc');
%Build & Solve
Opt = opti('f',f,'ineq',A,b,'bounds',lb,ub,'sos',sos,sosind,soswt,'options',opts)
[x,fval,exitflag,info] = solve(Opt)

%%
sos = [];
sos.type = sos; sos.index = sosind; sos.weight = soswt;
[x,fval,exitflag,info,Opt] = opti_mintprog(f,A,b,[],[],lb,ub,[],sos,opts)

%% BILP1
clc
%Objective & Constraints
f = -[6 5]';
A = [-3,5; 6,4; 3, -5; -6, -4]; 
b = [6;9;1;3];  
%Setup Options
opts = optiset('solver','glpk');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'int','BB','options',opts)
[x,fval,exitflag,info] = solve(Opt)
plot(Opt,3)

%%
[x,fval,exitflag,info,Opt] = opti_bintprog(f,A,b,[],[],[],opts)

%% BILP2
clc
%Objective & Constraints
f = -[9 5 6 4]';
A = [6 3 5 2; 0 0 1 1; -1 0 1 0; 0 -1 0 1];
b = [9; 1; 0; 0];
%Setup Options
opts = optiset('solver','glpk');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'int','BBBB','options',opts)
[x,fval,exitflag,info] = solve(Opt)

%%
[x,fval,exitflag,info,Opt] = opti_bintprog(f,A,b,[],[],[],opts)

%% QP1
clc
%Objective & Constraints
H = eye(3);
f = -[2 3 1]';
A = [1 1 1;3 -2 -3; 1 -3 2]; 
b = [1;1;1];
%Setup Options
opts = optiset('solver','ooqp');
%Build & Solve
Opt = opti('hess',H,'grad',f,'ineq',A,b,'options',opts)
[x,fval,exitflag,info] = solve(Opt)

%%
[x,fval,exitflag,info,lambda,Opt] = opti_quadprog(H,f,A,b,[],[],[],[],[],opts)

%% QP2
clc
%Objective & Constraints
H = [1 -1; -1 2];
f = -[2 6]';
A = [1 1; -1 2; 2 1];
b = [2; 2; 3];
%Setup Options
opts = optiset('solver','ooqp');
%Build & Solve
Opt = opti('qp',H,f,'ineq',A,b,'lb',[0;0],'options',opts)
[x,fval,exitflag,info] = solve(Opt)
%Plot
plot(Opt)

%%
[x,fval,exitflag,info,lambda,Opt] = opti_quadprog(H,f,A,b,[],[],[0;0],[],[],opts)

%% QP3
clc
%Objective & Constraints
H = [1 -1; -1 2];
f = -[2 6]';
A = [1 1; -1 2; 2 1];
b = [2; 2; 3];
Aeq = [1 1.5];
beq = 2;
%Setup Options
opts = optiset('solver','ooqp');
%Build & Solve
Opt = opti('qp',H,f,'eq',Aeq,beq,'ineq',A,b,'bounds',[0;0],[10;10],'options',opts)
[x,fval,exitflag,info] = solve(Opt)
%Plot
plot(Opt)

%%
[x,fval,exitflag,info,lambda,Opt] = opti_quadprog(H,f,A,b,Aeq,beq,[0;0],[10;10],[],opts)

%% SNLE1
clc
%Equations
f = @(x) [2*x(1) - x(2) - exp(-x(1));
          -x(1) + 2*x(2) - exp(-x(2))];
%Setup Options
x0 = [-5;5];
Opt = opti('fun',f,'x0',x0)
%Build & Solve
[x,fval,ef,stat] = solve(Opt)

%%
[x,fval,exitflag,info,Opt] = opti_fsolve(f,x0)

%% NLS1
clc
%Function
fun = @(x,xdata) x(1)*exp(x(2)*xdata);
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5];
%Setup Options
opts = optiset('solver','lmder','display','iter');
%Build & Solve
x0 = [100; -1]; % Starting guess
Opt = opti('fun',fun,'data',xdata,ydata,'ndec',2,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)
%Plot
plot(Opt,0.05,1)

%%
[x,fval,exitflag,info,Opt] = opti_lsqcurvefit(fun,x0,xdata,ydata,[],[],opts)

%% NLS2
clc
%Function
i = (1:40)';
fun = @(x) x(1)*exp(-x(2)*i) + x(3);
%Fitting Data
ydata=[5.8728, 5.4948, 5.0081, 4.5929, 4.3574, 4.1198, 3.6843, 3.3642, 2.9742, 3.0237, 2.7002, 2.8781,...
       2.5144, 2.4432, 2.2894, 2.0938, 1.9265, 2.1271, 1.8387, 1.7791, 1.6686, 1.6232, 1.571, 1.6057,...
       1.3825, 1.5087, 1.3624, 1.4206, 1.2097, 1.3129, 1.131, 1.306, 1.2008, 1.3469, 1.1837, 1.2102,...
       0.96518, 1.2129, 1.2003, 1.0743];
%Setup Options
opts = optiset('solver','nl2sol');
%Build & Solve
x0=[1.0; 0.0; 0.0];
Opt = opti('fun',fun,'ydata',ydata,'ndec',3,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)

%%
[x,fval,exitflag,info,Opt] = opti_lsqcurvefit(fun,x0,[],ydata,[],[],opts)

%% NLS3
clc
%Function
i = (1:40)';
fun = @(x) x(1)*exp(-x(2)*i) + x(3);
lb = [5.5;0;0];
%Fitting Data
ydata=[5.8728, 5.4948, 5.0081, 4.5929, 4.3574, 4.1198, 3.6843, 3.3642, 2.9742, 3.0237, 2.7002, 2.8781,...
       2.5144, 2.4432, 2.2894, 2.0938, 1.9265, 2.1271, 1.8387, 1.7791, 1.6686, 1.6232, 1.571, 1.6057,...
       1.3825, 1.5087, 1.3624, 1.4206, 1.2097, 1.3129, 1.131, 1.306, 1.2008, 1.3469, 1.1837, 1.2102,...
       0.96518, 1.2129, 1.2003, 1.0743];
%Setup Options
opts = optiset('solver','levmar','maxiter',5000);
%Build & Solve
x0=[1.0; 0.0; 0.0];
Opt = opti('fun',fun,'ydata',ydata,'lb',lb,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)

%%
[x,fval,exitflag,info,Opt] = opti_lsqcurvefit(fun,x0,[],ydata,lb,[],opts)

%% UNO Rosenbrock
clc
%Objective & Gradient
obj = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
%Setup Options
opts = optiset('solver','ipopt');
%Build & Solve
x0 = [0 0]';
Opt = opti('obj',obj,'ndec',2,'options',opts)
[x,fval,exitflag,info]= solve(Opt,x0)
%Plot
plot(Opt,[],1)

%%
[x,fval,exitflag,info,Opt] = opti_fminunc(obj,x0,opts)

%% NLP1 Hock & Schittkowski #71
clc
%Objective & Gradient
obj = @(x) x(1)*x(4)*sum(x(1:3)) + x(3);
grad = @(x) [ x(1)*x(4) + x(4)*sum(x(1:3));
              x(1)*x(4);
              x(1)*x(4) + 1;
              x(1)*sum(x(1:3)) ];          
%Linear Constraints
lb = ones(4,1);
ub = 5*ones(4,1);
%Nonlinear Constraints
nlcon = @(x) [ prod(x);
               sum(x.^2)];
nljac = @(x) [ prod(x)./x;
                2*x ];          
nlrhs = [25 40]';
nle = [1 0]'; % (>=, ==)
%Setup Options
opts = optiset('solver','ipopt','warnings','on','display','iter');
%Build & Solve
Opt = opti('obj',obj,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'options',opts)
x0 = [1 5 5 1]';
[x,fval,exitflag,info]= solve(Opt,x0)

%%
mwf = @(x) deal(obj(x),grad(x));
mwc = @(x) deal(-prod(x)+25,sum(x.^2)-40);

[x,fval,exitflag,info,Opt] = opti_fmincon(mwf,x0,[],[],[],[],lb,ub,mwc,opts)

%% NLP2 Hock & Schittkowski #38
clc
%Objective & Gradient
obj = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2 + 90*(x(4)-x(3)^2)^2 + (1-x(3))^2 + ...
      10.1*(x(2)-1)^2 + 10.1*(x(4)-1)^2 + 19.8*(x(2)-1)*(x(4)-1);  
grad = @(x) [ -400*x(1)*(x(2)-x(1)^2) - 2*(1-x(1));
              200*(x(2)-x(1)^2) + 20.2*(x(2)-1) + 19.8*(x(4)-1);
              -360*x(3)*(x(4)-x(3)^2) - 2*(1-x(3));
              180*(x(4)-x(3)^2) + 20.2*(x(4)-1) + 19.8*(x(2)-1)];
%Hessian (must not be tril for matlab!) & Hessian Structure       
hess = @(x) sparse(  [ 1200*x(1)^2-400*x(2)+2  -400*x(1)       0                          0
                       -400*x(1)               220.2           0                          19.8
                        0                      0               1080*x(3)^2- 360*x(4) + 2  -360*x(3)
                        0                      19.8            -360*x(3)                   200.2 ]);  
Hstr = @() sparse([ 1  0  0  0 
                    1  1  0  0
                    0  0  1  0
                    0  1  1  1 ]);                               
%Constraints
lb = [-10 -10 -10 -10]';
ub = [10 10 10 10]';
%Setup Options
opts = optiset('solver','ipopt');
%Build & Solve
Opt = opti('obj',obj,'grad',grad,'hess',hess,'Hstr',Hstr,'bounds',lb,ub,'options',opts)
x0 = [-3  -1  -3  -1]';
[x,fval,exitflag,info]= solve(Opt,x0)

%%
mwf = @(x) deal(obj(x),grad(x));
mwc = [];

[x,fval,exitflag,info,Opt] = opti_fmincon(mwf,x0,[],[],[],[],lb,ub,mwc,opts)

%% NLP3 Hock & Schittkowski #51
clc
%Objective & Gradient
obj = @(x) (x(1) - x(2))^2 + (x(2) + x(3) - 2)^2 + (x(4) - 1)^2 + (x(5) - 1)^2;
grad = @(x) 2*[ x(1) - x(2);
                x(2) + x(3) - 2 - x(1) + x(2);
                x(2) + x(3) - 2;
                x(4) - 1;
                x(5) - 1 ];
%Nonlinear Constraints & Jacobian Structure
nlcon = @(x) [ x(1) + 3*x(2);
               x(3) + x(4) - 2*x(5);
               x(2) - x(5) ];
nljac = @(x) sparse([ 1  3  0  0  0;
	                  0  0  1  1 -2;
	                  0  1  0  0 -1 ]);
nljacstr = @() sparse([1 1 0 0 0;
                       0 0 1 1 1;
                       0 1 0 0 1]);
nlrhs = [4 0 0]';
nle = [0 0 0]';
%Setup Options
opts = optiset('solver','nlopt');
%Build & Solve
Opt = opti('obj',obj,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'nljacstr',nljacstr,'options',opts)
x0 = [ 2.5 0.5 2 -1 0.5 ];
[x,fval,exitflag,info]= solve(Opt,x0)

%%
mwf = @(x) deal(obj(x),grad(x));
mwc = @(x) deal([],[x(1) + 3*x(2)-4;
                    x(3) + x(4) - 2*x(5);
                    x(2) - x(5)]);

[x,fval,exitflag,info,Opt] = opti_fmincon(mwf,x0,[],[],[],[],[],[],mwc,opts)
