%% Test Problems for NOMAD Interface
% J.Currie 2012

% Many problems are from Nick Sahinidis's BARON Global Test Set and sample
% solutions are obtained via BARON.

%% Options File
%Options which are available via MEX interface. The remainder are available
%via a parameters file, supplied via 'paramfile' in nomadset.

clc
nomadset

%% NOMAD Info
clc
nomad('-i') %info
nomad('-v') %ver

%% Rosenbrock [x = 1,1, fval = 0]
clc
fun = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
x0 = [0 0]';

[x,fval,ef,iter] = nomad(fun,x0)

%% Constrained Rosenbrock [x = .9488,.9, fval = .0026]
clc
fun = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
x0 = [0 0]';
lb = [0;0];
ub = [1;0.9];

[x,fval,ef,iter] = nomad(fun,x0,lb,ub)

%% Rosenbrock with Inf Bound [x = 1,1, fval = 0]
clc
fun = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
x0 = [-2 1]';
lb = [-1;-1.5];
ub = [100;100];
opts = nomadset('display_degree',2,'min_mesh_size','1e-007','max_time',5);

[x,fval,ef,iter] = nomad(fun,x0,lb,ub,[],[],'cc',opts)

%% St_e01 [x = 6,0.6667, fval = -6.6667]
clc
%Objective
fun = @(x) -x(1) - x(2);
%Nonlinear Constraints
nlcon = @(x) x(1)*x(2);
nlrhs = 4;
%Bounds
lb = [0;0];
ub = [6;4];
x0 = [1,1]';
opts = nomadset('display_degree',2,'min_mesh_size','1e-007','direction_type','ortho 2n');
tic
[x,fval,ef,iter,nfval] = nomad(fun,x0,lb,ub,nlcon,nlrhs,'CC',opts)
toc
%% Basic [x = 3.5355, 3.5355, fval = 12.8579]
clc
%Objective 
fun = @(x) sum((x-1).^2);
%Constraints
nlcon = @(x) [sum((x-1).^2)-25; 25-sum(x.^2)];
nlrhs = [0;0];
lb = [-10;-10];
ub = [10;10];
x0 = [0;0];

[x,fval,ef,iter] = nomad(fun,x0,lb,ub,nlcon,nlrhs,[],opts)

%% 4x Rosenbrock [1,1,1,1]
clc
%Objective
fun = @(x) sum(100*(x(1:end-1).^2-x(2:end)).^2+(x(1:end-1)-1).^2);
x0 = randn(4,1);
[x,fval,ef,iter] = nomad(fun,x0,[],[],[],[],[],opts)

%% Canoe [x = -0.014, 0.1241, fval = 849.1080) (needs problem data)
% clc
% fun = @(X) norm(ss(A+B*diag(X)*C, B1+B*diag(X)*D21,C1+D12*diag(X)*C,D11+D12*diag(X)*D21),inf,10^-7)*(max(real(eig(A+B*diag(X)*C)) < 0)) + (max(real(eig(A+B*diag(X)*C)))>=0)*10^100;
% 
% x0 = [-0.01274; 0.12034];
% lb=[-0.015;0.12];
% ub=[-0.010;0.13];
% 
% [x,fval,ef,iter] = nomad(fun,x0,lb,ub,[],[],[],opts)

%% St_e08 [x = 0.1294, 0.4830, fval = 0.7418]
clc
%Objective
fun = @(x) 2*x(1)+x(2);
%Constraints
nlcon = @(x) [-16*x(1)*x(2); (-4*x(1)^2) - 4*x(2)^2];
nlrhs = [-1;-1];
lb = [0;0];
ub = [1;1];
x0 = [0;0];

[x,fval,ef,iter] = nomad(fun,x0,lb,ub,nlcon,nlrhs,[],opts)

%% St_e09 [x = 0.5, 0.5, fval = -0.5]
clc
%Objective
fun = @(x) -2*x(1)*x(2);
%Constraints
nlcon = @(x) 4*x(1)*x(2) + 2*x(1) + 2*x(2);
nlrhs = 3;
lb = [0;0];
ub = [1;1];
x0 = [0.4;0.4];

[x,fval,ef,iter] = nomad(fun,x0,lb,ub,nlcon,nlrhs,[],opts)

%% St_e17 [x = 8.1651, 7.5685, fval = 376.2889]
clc
%Objective
fun = @(x) 29.4*x(1) + 18*x(2);
%Constraints
nlcon = @(x) -(x(1) - 0.2458*x(1)^2/x(2));
nlrhs = -6;
lb = [0;1e-5];
ub = [115.8;30];
x0 = [0;1e-5];

[x,fval,ef,iter] = nomad(fun,x0,lb,ub,nlcon,nlrhs,[],opts)

%% St_e18 [x = -1.4142, -1.4142, fval = -2.8284]
clc
%Objective
fun = @(x) x(1) + x(2);
%Constraints
nlcon = @(x) [-(x(1)^2 + x(2)^2); x(1)^2 + x(2)^2; -x(1) + x(2); x(1) - x(2)];
nlrhs = [-1;4;1;1];
lb = [-2;-2];
ub = [2;2];
x0 = [0;0];

[x,fval,ef,iter] = nomad(fun,x0,lb,ub,nlcon,nlrhs,[],opts)

%% St_e19 [x = -3.1736, 1.7245, fval = -118.7049]
clc
%Objective
fun = @(x) x(1)^4 - 14*x(1)^2 + 24*x(1) - x(2)^2;
%Constraints
nlcon = @(x) [-x(1) + x(2); (-x(1)^2) - 2*x(1) + x(2)];
nlrhs = [8;-2];
lb = [-8;0];
ub = [10;10];
x0 = [0;0];

[x,fval,ef,iter] = nomad(fun,x0,lb,ub,nlcon,nlrhs,[],opts)

%% ex2_1_3 [fval = -15] (takes a while)
clc
%Objective
fun = @(x) 5*x(1) - 0.5*(10*x(1)*x(1) + 10*x(2)*x(2) + 10*x(3)*x(3) + 10*x(4)*x(4)) + 5*x(2) + 5*x(3) + 5*x(4) - x(5) - x(6) - x(7) - x(8) - x(9) - x(10) - x(11) - x(12) - x(13);
%Constraints
nlcon = @(x) [2*x(1) + 2*x(2) + x(10) + x(11);
              2*x(1) + 2*x(3) + x(10) + x(12);
              2*x(2) + 2*x(3) + x(11) + x(12)
              -8*x(1) + x(10);
              -8*x(2) + x(11);
              -8*x(3) + x(12);
              -2*x(4) - x(5) + x(10);
              -2*x(6) - x(7) + x(11);
              -2*x(8) - x(9) + x(12)];
nlrhs = [10;10;10;0;0;0;0;0;0];
lb = zeros(13,1);
ub = 1*ones(13,1); ub(10:12) = 5;
x0 = lb;

[x,fval,ef,iter] = nomad(fun,x0,lb,ub,nlcon,nlrhs,[],opts)

%% Wolfram problem
clc
%Objective
obj = @(x) [x(1) - sin(2*x(1) + 3*x(2)) - cos(3*x(1) - 5*x(2));
          x(2) - sin(x(1) - 2*x(2)) + cos(x(1) + 3*x(2))];
fun = @(x) norm(obj(x));

lb = [-4;-4];
ub = [4;4];
x0 = [0;-3];

opts = nomadset('display_degree',2,'vns_search',0.9,'f_target',1e-6);

[xr,fval,ef,iter] = nomad(fun,x0,lb,ub,[],[],[],opts)


%% Test Results for various Solvers
% [xr2,fval2,ef,iter] = pswarm(fun,x0,lb,ub)
% [xr3,fval3,ef,iter] = ga(fun,2,[],[],[],[],lb,ub)
% [xr4,fval4,ef,iter] = patternsearch(fun,x0,[],[],[],[],lb,ub)
% [xr5,fval5,ef,iter] = simulannealbnd(fun,x0,lb,ub)
% [xr6,fval6,ef,iter] = opti_fmincon(fun,x0,[],[],[],[],lb,ub)
% 
% prob = optiprob('fun',fun,'bounds',lb,ub,'x0',x0);
% opts = optiset('solver','nlopt','solverOpts',nloptset('algorithm','GN_DIRECT'));
% Opt = opti(prob,opts);
% [xr7,fval7] = solve(Opt);
% 
% % n = 3e2; xmin = -4; xmax = 4;
% % x = linspace(xmin,xmax,n);
% % y = linspace(xmin,xmax,n);
% % Z = zeros(n,n);
% % for i = 1:n
% %     for j = 1:n
% %         Z(j,i) = fun([x(i),y(j)]);
% %     end
% % end
% surf(x,y,Z)
% colormap summer
% shading flat
% view(0,90);
% xlabel('x1'); ylabel('x2'); title('Wolfram Global Problem');      
% 
% hold on
% plot3(x0(1),x0(2),10,'kx','markersize',15);
% plot3(xr(1),xr(2),5,'ro'); text(xr(1)+0.1,xr(2),5,sprintf('nmd %4.4f',fval))
% plot3(xr2(1),xr2(2),5,'bo'); text(xr2(1)+0.1,xr2(2),5,sprintf('swm %4.4f',fval2))
% plot3(xr3(1),xr3(2),5,'mo'); text(xr3(1)+0.1,xr3(2)-0.2,5,sprintf('ga %4.4f',fval3))
% plot3(xr4(1),xr4(2),5,'co'); text(xr4(1)+0.1,xr4(2),5,sprintf('pat %4.4f',fval4));
% plot3(xr5(1),xr5(2),5,'ko'); text(xr5(1)+0.1,xr5(2),5,sprintf('anl %4.4f',fval5));
% plot3(xr6(1),xr6(2),5,'yo'); text(xr6(1)+0.1,xr6(2),5,sprintf('ipt %4.4f',fval6));
% plot3(xr7(1),xr7(2),5,'wo'); text(xr7(1)+0.1,xr7(2)+0.2,5,sprintf('npt %4.4f',fval7));
% hold off;

%% MINLP 1 [fval = -5]
clc
fun = @(x) (x(1) - 5)^2 + x(2)^2 - 25;
nlcon = @(x) x(1)^2 - x(2) + 0.5;
nlrhs = 0;     
xtype = 'II';
x0 = [0;0];

opts = nomadset('display_degree',2);

[xr,fval,ef,iter] = nomad(fun,x0,[],[],nlcon,nlrhs,xtype,opts)

%% MINLP2 [fval = -1.5]
clc
fun = @(x) -x(1) - x(2) - x(3);
nlcon = @(x) [ (x(2) - 1./2.)*(x(2) - 1./2.) + (x(3) - 1./2.)*(x(3) - 1./2.);
                x(1) - x(2);
                x(1) + x(3) + x(4)];          
nlrhs = [1/4;0;2];
ub = [1;10;10;5];
lb = [0;0;0;0];
xtype = 'BCCI';
x0 = [0;0;0;0];

[xr,fval,ef,iter] = nomad(fun,x0,lb,ub,nlcon,nlrhs,xtype,opts)

%% MINLP 3 [fval = 424.1355]
clc
fun = @(x) 20 + x(1)^2 + x(2)^2 - 10*(cos(2*pi*x(1)) + cos(2*pi*x(2)));
lb = [5*pi;-20*pi];
ub = [20*pi;-4*pi];
xtype = 'IC';
x0 = [30;-30]

[xr,fval,ef,iter] = nomad(fun,x0,lb,ub,[],[],xtype,opts)

%% BI Objective (NOMAD Example)
clc
fun = @(x) [x(5); sum((x-1).^2)-25];
nlcon = @(x) 25-sum((x+1).^2);          
nlrhs = 0;
lb = [-6;-6;-6;-6;-6];
ub = [5;6;7;Inf;Inf];
xtype = 'CCCCC';
x0 = [0;0;0;0;0];
opts = nomadset('display_degree',2,'multi_overall_bb_eval',100,'bb_output_type',{'EB'});

[xr,fval,ef,iter] = nomad(fun,x0,lb,ub,nlcon,nlrhs,xtype,opts)


%% Test Cell Based Options
clc
fun = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
x0 = [0 0]';

opts = nomadset('display_degree',2,'initial_mesh_size',{'0 10','1 10'});

[xr,fval,ef,iter] = nomad(fun,x0,[],[],[],[],[],opts)


%% Test Cell based Input
clc
fun = @(x) (x(1) - 5)^2 + x(2)^2 - 25;
nlcon = @(x) x(1)^2 - x(2) + 0.5;
            
x0 = [0;0];

opts = nomadset('display_degree',2,'bb_input_type',{'0 I','1 I'});

[xr,fval,ef,iter] = nomad(fun,x0,[],[],nlcon,[],[],opts)