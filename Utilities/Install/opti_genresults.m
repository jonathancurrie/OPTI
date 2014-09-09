%% Generate Test Problems + Results for opti_Install_Test
% Assumes you are in base opti directory!
% All simple, quick problems

%% Linear Programming
clc
clear
no = 5;
lp_tprob = cell(no,1);
lp_sval = zeros(no,1);
i = 1;

%LP 1
f = -[-1, 2]';
A = [2, 1;-4, 4];
b = [5, 5]';  
lp_tprob{i} = optiprob('f',f,'ineq',A,b); 
lp_sval(i) = -3.75;
i = i + 1;

%LP 2
f = -[50, 100];
A = [10, 5;4, 10; 1, 1.5];
b = [2500, 2000, 450]';
lp_tprob{i} = optiprob('f',f,'ineq',A,b);            
lp_sval(i) = -21875;
i = i + 1;
        
%LP 3
f = [2, 3, 7, 7];
A = [1, 1, -2, -5;-1, 2, 1, 4];
b = [2, -3]';
e = [1, 1];
lb = zeros(4,1); 
ub = [30 100 20 1]';   
lp_tprob{i} = optiprob('f',f,'mix',A,b,e,'bounds',lb,ub);            
lp_sval(i) = 4;
i = i + 1;

%LP 4
f = [1, 2, 3, 7, 8, 8];
A = [5, -3, 2, -3, -1, 2; -1, 0, 2, 1, 3, -3;1, 2, -1, 0, 5, -1];
b = [-5, -1, 3]';
e = [1, 1, 1];
lb = zeros(6,1);
ub = 10*ones(6,1);   
lp_tprob{i} = optiprob('f',f,'mix',A,b,e,'bounds',lb,ub);            
lp_sval(i) = 3;
i = i + 1;
        
%LP 5 
f = -[1 2 3]';
A = [-1,1,1; 1,-3,1];
b = [20,30]';
Aeq = [1 1 1];
beq = 40;
lb = [0 0 0]';
ub = [40 inf inf]';
lp_tprob{i} = optiprob('grad',f,'ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub);
lp_sval(i) = -97.5;
i = i + 1;

%LP Clean Up + Save
clear f A b e Aeq beq lb ub i no
save 'Utilities/Install/Test Results/lp_test_results.mat'


%% Mixed Integer Linear Programming
clc
clear
no = 5;
milp_tprob = cell(no,1);
milp_sval = zeros(no,1);
i = 1;

%MILP 1
f = -[-1, 2]';
A = [2, 1;-4, 4];
b = [5, 5]';  
xint = 'II';
milp_tprob{i} = optiprob('f',f,'ineq',A,b,'int',xint); 
milp_sval(i) = -3;
i = i + 1;

%MILP 2
f = [2, 3, 7, 7];
A = [1, 1, -2, -5;-1, 2, 1, 4];
b = [2, 3]';
e = [1, 1];
lb = zeros(4,1); 
ub = [30 100 20 1]';  
xint = 'CICI';
milp_tprob{i} = optiprob('f',f,'mix',A,b,e,'bounds',lb,ub,'int',xint);
milp_sval(i) = 6;
i = i + 1;

%MILP 3
f = -[592, 381, 273, 55, 48, 37, 23];
A = [3534, 2356, 1767, 589, 528, 451, 304];
b = 119567;
lb = zeros(7,1); 
ub = [100 50 33 20 77 44 20]';
xint = 'IIIIIII';
milp_tprob{i} = optiprob('f',f,'ineq',A,b,'bounds',lb,ub,'int',xint);            
milp_sval(i) = -19979;
i = i + 1;

%MILP 4
f = -[6 5]';
A = [-3,5; 6,4; 3, -5; -6, -4]; 
b = [6;9;1;3];  
xint = 'BC';
milp_tprob{i} = optiprob('grad',f,'ineq',A,b,'int',xint);
milp_sval(i) = -9.75;
i = i + 1;

%MILP 5
f = -[1 2 3 1]'; 
A = [-1 1 1 10; 1 -3 1 0]; 
b = [20;30];  
Aeq = [0 1 0 -3.5];
beq = 0;
lb = [0 0 0 2]';
ub = [40 inf inf 3]';
xint = 'CCCI';
%Build & Options
milp_tprob{i} = optiprob('grad',f,'ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub,'int',xint);
milp_sval(i) = -122.5;

clear f A b e Aeq beq lb ub xint i no
save 'Utilities/Install/Test Results/milp_test_results.mat'

%% Quadratic Programming
clc
clear
no = 3;
qp_tprob = cell(no,1);
qp_sval = zeros(no,1);
i = 1;

%QP 1
H = eye(3);
f = -[2 3 1]';
A = [1 1 1;3 -2 -3; 1 -3 2]; 
b = [1;1;1];  
qp_tprob{i} = optiprob('H',H,'f',f,'ineq',A,b); 
qp_sval(i) = -2.83333333333333;
i = i + 1;

%QP 2
H = [1 -1; -1 2];
f = -[2 6]';
A = [1 1; -1 2; 2 1];
b = [2; 2; 3];
lb = [0;0];
qp_tprob{i} = optiprob('qp',H,f,'ineq',A,b,'lb',lb);
qp_sval(i) = -8.22222220552525;
i = i + 1;

%QP 3
H = [1 -1; -1 2];
f = -[2 6]';
A = [1 1; -1 2; 2 1];
b = [2; 2; 3];
Aeq = [1 1.5];
beq = 2;
lb = [0;0];
ub = [10;10];
qp_tprob{i} = optiprob('qp',H,f,'eq',Aeq,beq,'ineq',A,b,'bounds',lb,ub);
qp_sval(i) = -6.41379310344827;
i = i + 1;

clear H f A b e Aeq beq lb ub i no
save 'Utilities/Install/Test Results/qp_test_results.mat'

%% MI Quadratic Programming
clc
clear
no = 2;
miqp_tprob = cell(no,1);
miqp_sval = zeros(no,1);
i = 1;

%MIQP 1
H = [1 -1; -1 2];
f = -[2 6]';
A = [1 1; -1 2];
b = [3; 3.5];
miqp_tprob{i} = optiprob('hess',H,'grad',f,'ineq',A,b,'bounds',[0;0],[10;10],'int',1);
miqp_sval(i) = -11.5;
i = i + 1;

%MIQP 2
H = eye(3);
f = -[2 3 1]';
A = [1 1 1;3 -2 -3; 1 -3 2]; 
b = [1;1;1];
miqp_tprob{i} = optiprob('hess',H,'grad',f,'ineq',A,b,'int','CIC');
miqp_sval(i) = -2.75;
i = i + 1;

clear H f A b e Aeq beq lb ub i no
save 'Utilities/Install/Test Results/miqp_test_results.mat'

%% SDP Programming
clc
clear
no = 5;
sdp_tprob = cell(no,1);
sdp_sval = zeros(no,1);
i = 1;

%SDP1
f = 1;
A = eye(2);
C = -[0 sqrt(2); sqrt(2) 0];
sdp = sparse([C(:) A(:)]);
sdp_tprob{i} = optiprob('f',f,'sdcone',sdp);
sdp_sval(i) = sqrt(2);
i = i + 1;
sdp = [];

%SDP2
f = [1 1]';
C = -[0 2; 2 0];
A0 = [1 0; 0 0];
A1 = [0 0; 0 1];
sdp = sparse([C(:) A0(:) A1(:)]);
lb = [0;0];
ub = [10;10];
sdp_tprob{i} = optiprob('f',f,'sdcone',sdp,'bounds',lb,ub);             
sdp_sval(i) = 4;
i = i + 1;
sdp = [];

%SDP3
f = [1 0 0 0]';
%SDP Constraint1 [x2 x3;x3 x4] <= x1*eye(2)
C = zeros(2);
A0 = eye(2);
A1 = -[1 0; 0 0];
A2 = -[0 1; 1 0];
A3 = -[0 0; 0 1];
sdp{1} = sparse([C(:) A0(:) A1(:) A2(:) A3(:)]);
%SDP Constraint2 [x2 x3; x3 x4] >= [1 0.2; 0.2 1]
C = [1 0.2; 0.2 1];
A0 = zeros(2);
A1 = [1 0; 0 0];
A2 = [0 1; 1 0];
A3 = [0 0; 0 1];
sdp{2} = sparse([C(:) A0(:) A1(:) A2(:) A3(:)]);
sdp_tprob{i} = optiprob('f',f,'sdcone',sdp);          
sdp_sval(i) = 1.2;   
i = i + 1;
sdp = [];

%SDP4
f = [1 1 1]';
C = -[0 1 2; 1 0 3; 2 3 100];
A0 = [1 0 0; 0 0 0; 0 0 0];
A1 = [0 0 0; 0 1 0; 0 0 0];
A2 = zeros(3);
sdp = sparse([C(:) A0(:) A1(:) A2(:)]);
lb = [10;0;0];
ub = [1000;1000;1000];
sdp_tprob{i} = optiprob('f',f,'sdcone',sdp,'bounds',lb,ub);             
sdp_sval(i) = 10.1787148490962;
i = i + 1;
sdp = [];

%SDP5
f = -[1 2];
%SDP Constraint 1
C = [2 1; 1 2];
A0 = -[3 1; 1 3];
A1 = -zeros(2);
sdp{1} = sparse([C(:) A0(:) A1(:)]);
%SDP Constraint 2
C = [3 0 1; 0 2 0; 1 0 3];
A0 = -zeros(3);
A1 = -[3 0 1; 0 4 0; 1 0 5];
sdp{2} = sparse([C(:) A0(:) A1(:)]);
%SDP Constraint 3
C = zeros(2);
A0 = -[1 0; 0 0];
A1 = -[0 0; 0 1];
sdp{3} = sparse([C(:) A0(:) A1(:)]);
sdp_tprob{i} = optiprob('f',f,'sdcone',sdp);          
sdp_sval(i) = 2.74999999779937;

clear f A0 A1 A2 A3 sdp lb ub i no
save 'Utilities/Install/Test Results/sdp_test_results.mat'

%% Nonlinear Least Squares
clc
clear
no = 5;
nls_tprob = cell(no,1);
nls_sval = zeros(no,1);
i = 1;

%NLS 1
fun = @(x) [10*(x(2)-x(1)^2); 1 - x(1)];
ydata = zeros(2,1);
x0 = [-1.2;1];
nls_tprob{i} = optiprob('fun',fun,'ydata',ydata,'x0',x0);            
nls_sval(i) = 0;
i = i + 1;

%NLS 2
fun = @(x) [1.5 - x(1)*(1 - x(2)); 2.25 - x(1)*(1 - x(2)^2); 2.625 - x(1)*(1 - x(2)^3) ];
ydata = zeros(3,1);
x0 = [1;1];
nls_tprob{i} = optiprob('fun',fun,'ydata',ydata,'x0',x0);            
nls_sval(i) = 0;
i = i + 1;

%NLS 3
fun = @(x) [1e4*x(1)*x(2) - 1;
            exp(-x(1)) + exp(-x(2)) - 1.0001];
ydata = zeros(2,1);
x0 = [0;1];
nls_tprob{i} = optiprob('fun',fun,'ydata',ydata,'x0',x0);            
nls_sval(i) = 0;
i = i + 1;

%NLS 4
ii = (1:10)';
fun = @(x) 2 + 2*ii - (exp(0.2578*ii) + exp(0.2578*ii));
ydata = zeros(10,1);
x0 = [0.3;0.4];
nls_tprob{i} = optiprob('fun',fun,'ydata',ydata,'x0',x0);            
nls_sval(i) = 124.36226865;
i = i + 1;

%NLS 5
fun = @(x) [10*(x(2)-x(1)^2);
            1 - x(1);
            sqrt(90)*(x(4)-x(3)^2);
            1 - x(3);
            sqrt(10)*(x(2) + x(4) - 2);
            10^(-0.5)*(x(2)-x(4))]; 
x0 = [3;-1;-3;-1];
ydata = zeros(6,1);
nls_tprob{i} = optiprob('fun',fun,'ydata',ydata,'x0',x0);            
nls_sval(i) = 0;

clear fun ydata x0 ii i no
save 'Utilities/Install/Test Results/nls_test_results.mat'

%% Nonlinear Programming
clc
clear
no = 3;
nlp_tprob = cell(no,1);
nlp_sval = zeros(no,1);
i = 1;

%NLP 1
fun = @(x) 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
grad = @(x)[[2*x(1)-400*x(1)*(-x(1)^2+x(2))-2];[-200*x(1)^2+200*x(2)]];
lb = [-inf; -1.5];
ub = [100; 100];
x0 = [-2; 1];        
nlp_tprob{i} = optiprob('obj',fun,'grad',grad,'bounds',lb,ub,'x0',x0);            
nlp_sval(i) = 0;
i = i + 1;

%NLP 2
fun = @(x) 1/3*(x(1) + 1)^3 + x(2);
grad = @(x)[[(x(1)+1)^2];[1]];
lb = [1;0];
ub = [100;100];
x0 = [1.125;0.125];
nlp_tprob{i} = optiprob('obj',fun,'grad',grad,'bounds',lb,ub,'x0',x0);
nlp_sval(i) = 8/3;
i = i + 1;

%NLP 3
fun = @(x) sin(x(1) + x(2)) + (x(1) - x(2))^2 - 1.5*x(1) + 2.5*x(2) + 1; 
grad = @(x)[[2*x(1)-2*x(2)+cos(x(1)+x(2))-1.5];[2*x(2)-2*x(1)+cos(x(1)+x(2))+2.5]];
lb = [-1.5;-3];
ub = [4;3];
x0 = [0;0];
nlp_tprob{i} = optiprob('obj',fun,'grad',grad,'bounds',lb,ub,'x0',x0);
nlp_sval(i) = -0.5*sqrt(3)-pi/3;
i = i + 1;

clear fun grad lb ub x0 i no
save 'Utilities/Install/Test Results/nlp_test_results.mat'

%% MI Nonlinear Programming
clc
clear
no = 1;
minlp_tprob = cell(no,1);
minlp_sval = zeros(no,1);
i = 1;

%MINLP 1
obj = @(x) (x(1) - 5)^2 + x(2)^2 - 25;
nlcon = @(x) -x(1)^2 + x(2)-0.5;
nlrhs = 0;
nle = 1;  
x0 = [4;0];
minlp_tprob{i} = optiprob('obj',obj,'ndec',2,'ivars',[1 2],'nlmix',nlcon,nlrhs,nle,'x0',x0);
minlp_sval(i) = -5;
i = i + 1;

clear obj A b nlcon nlrhs nle lb ub x0 i no
save 'Utilities/Install/Test Results/minlp_test_results.mat'