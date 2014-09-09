%% MINLP and NLP Testing
% A collection of NLPs from Hock-Schittkowski and other sources with
% solutions for testing the SCIP - MATLAB Interface.

% Run the entire file by using F5

%Check SCIP is available
if(~checkSolver('SCIP',0))
    fprintf(2,'\n\nSCIP IS NOT AVAILABLE, NO SCIP TESTS CAN BE RUN\n\n');
    return;
end


%%

i = 1;
clear fval fmin
opts = optiset('display','iter');

%% NLP1
clc
a = 100;
fun = @(x) a*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
lb = [-inf; -1.5];
x0 = [-2; 1];                   
fmin(i) = 0;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,[],[],[],[],[],x0,opts);
i = i + 1;

%% NLP2
fun = @(x) 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
lb = [-inf; 1.5];
x0 = [-2; 1];        
fmin(i) = 0.0504261879;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,[],[],[],[],[],x0,opts);
i = i + 1;

%% NLP3
fun = @(x) x(2) + 1e-5*(x(2)-x(1))^2;
lb = [-inf;0];
x0 = [10;1];
fmin(i) = 0;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,[],[],[],[],[],x0,opts);
i = i + 1;

%% NLP4
fun = @(x) 1/3*(x(1) + 1)^3 + x(2);
lb = [1;0];
x0 = [1.125;0.125];
fmin(i) = 8/3;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,[],[],[],[],[],x0,opts);
i = i + 1;

%% NLP5 (won't run - trig)
% fun = @(x) sin(x(1) + x(2)) + (x(1) - x(2))^2 - 1.5*x(1) + 2.5*x(2) + 1;    
% lb = [-1.5;-3];
% ub = [4;3];
% x0 = [0;0];
% fmin(i) = -0.5*sqrt(3)-pi/3;
% 
% [~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,ub,[],[],[],[],x0,opts);
% i = i + 1;

%% NLP6
fun = @(x) (1-x(1))^2;
nlcon = @(x) 10*(x(2) - x(1)^2);
cl = 0;
cu = 0;
x0 = [-1.2;1];
fmin(i) = 0;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP7
fun = @(x) log(1+x(1)^2) - x(2);
nlcon = @(x) (1 + x(1)^2)^2 + x(2)^2 - 4;
cl = 0;
cu = 0;
x0 = [2;2];
fmin(i) = -sqrt(3);

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP8 
clc
fun = @(x) -1;
nlcon = @(x) [x(1)^2 + x(2)^2 - 25;
              x(1)*x(2)-9];
cl = [0;0];
cu = [0;0];
x0 = [2;1];
fmin(i) = -1;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP9 (won't run - trig)
% fun = @(x) sin(pi*x(1)/12)*cos(pi*x(2)/16);
% nlcon = @(x) 4*x(1)-3*x(2);
% cl = 0;
% cu = 0;
% x0 = [0;0];
% fmin(i) = -.5;
% 
% [~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);
% i = i + 1;

%% NLP10
fun = @(x) x(1) - x(2);
nlcon = @(x) -3*x(1)^2 + 2*x(1)*x(2) - x(2)^2 + 1;
cl = 0;
cu = inf;
x0 = [-10;10];
fmin(i) = -1;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP11 
fun = @(x) (x(1) - 5)^2 + x(2)^2 - 25;
nlcon = @(x) -x(1)^2 + x(2);
cl = 0;
cu = inf;
x0 = [4.9;0.1];
fmin(i) = -8.498464223;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP12
fun = @(x) 0.5*x(1)^2 + x(2)^2 - x(1)*x(2) - 7*x(1) - 7*x(2);
nlcon = @(x) 25-4*x(1)^2-x(2)^2;
cl = 0;
cu = inf;
x0 = [0;0];
fmin(i) = -30;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP13
fun = @(x) (x(1)-2)^2 + x(2)^2;
nlcon = @(x)(1-x(1))^3 - x(2);
cl = 0;
cu = inf;
lb = [0;0];
x0 = [-2;-2];
fmin(i) = 0.998831;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP14
fun = @(x) (x(1)-2)^2 + (x(2)-1)^2;
nlcon = @(x) [-0.25*x(1)^2 - x(2)^2 + 1;
               x(1) - 2*x(2) + 1];
cl = [0;0];
cu = [inf;0];
x0 = [2;2];
fmin(i) = 9-2.875*sqrt(7);

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP15
fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
nlcon = @(x) [x(1)*x(2)-1
              x(1) + x(2)^2];
cl = [0;0];
cu = [inf;inf];
ub = [0.5;inf];
x0 = [-2,1];
fmin(i) = 306.5;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],ub,nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP16
fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
nlcon = @(x) [x(1) + x(2)^2
              x(1)^2 + x(2)];
cl = [0;0];
cu = [inf;inf];
lb = [-0.5;-inf];
ub = [0.5;1];
x0 = [-2,1];
fmin(i) = 0.25;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,ub,nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP17
fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
nlcon = @(x) [x(2)^2 - x(1)
              x(1)^2 - x(2)];
cl = [0;0];
cu = [inf;inf];
lb = [-0.5;-inf];
ub = [0.5;1];
x0 = [-2,1];
fmin(i) = 1;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,ub,nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP18
fun = @(x) 0.01*x(1)^2 + x(2)^2;
nlcon = @(x) [x(1)*x(2) - 25;
              x(1)^2 + x(2)^2 - 25];
cl = [0;0];
cu = [inf;inf];
lb = [2;0];
ub = [50;50];
x0 = [2;2];
fmin(i) = 5;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,ub,nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP 19 %PQR-T1-4 Betts[8], Gould[27]
fun = @(x) (x(1)-10)^3 + (x(2)-20)^3;
nlcon = @(x) [(x(1)-5)^2 + (x(2)-5)^2 - 100;
               -(x(2)-5)^2 - (x(1)-6)^2 + 82.81];
cl = [0;0];
cu = [inf;inf];
lb = [13;0];
ub = [100;100];
x0 = [20.1;5.84];
fmin(i) = -6961.81381;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,ub,nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP20
fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
nlcon = @(x) [x(1) + x(2)^2;
              x(1)^2 + x(2);
              x(1)^2 + x(2)^2 - 1];
cl = [0;0;0];
cu = [inf;inf;inf];
lb = [-0.5;-inf];
ub = [0.5;inf];
x0 = [-2;1];
fmin(i) = 81.5-25*sqrt(3);

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,ub,nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP21
fun = @(x) 0.01*x(1)^2 + x(2)^2 - 100;
nlcon = @(x) 10*x(1) - x(2) - 10; %note linear
cl = 0;
cu = inf;
lb = [2;-50];
ub = [50;50];
x0 = [-1;-1];
fmin(i) = -99.96;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,ub,nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP22
fun = @(x) (x(1)-2)^2 + (x(2)-1)^2;
nlcon = @(x) [-x(1) - x(2) + 2; %note linear
              -x(1)^2 + x(2)];
cl = [0;0];
cu = [inf;inf];
x0 = [2;2];
fmin(i) = 1;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP23
fun = @(x) x(1)^2 + x(2)^2;
nlcon = @(x) [x(1) + x(2) - 1;   
              x(1)^2 + x(2)^2 - 1;
              9*x(1)^2 + x(2)^2 - 9;
              x(1)^2 - x(2);
              x(2)^2 - x(1)];
cl = [0;0;0;0;0];
cu = [inf;inf;inf;inf;inf];
lb = [-50;-50];
ub = [50;50];
x0 = [3;1];
fmin(i) = 2;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,ub,nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP24
fun = @(x) 1/(27*sqrt(3)) * ((x(1) - 3)^2 - 9)*x(2)^3;
nlcon = @(x) [x(1)/sqrt(3) - x(2);
              x(1) + sqrt(3)*x(2);
              -x(1) - sqrt(3)*x(2) + 6];
cl = [0;0;0];
cu = [inf;inf;inf];
lb = [0;0];
x0 = [1;0.5];
fmin(i) = -1;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP25 (not in anonymous function form)

%% NLP26
fun = @(x) (x(1)-x(2))^2 + (x(2)-x(3))^4;
nlcon = @(x) (1+x(2)^2)*x(1) + x(3)^4 - 3;
cl = 0;
cu = 0;
x0 = [-2.6;2;2];
fmin(i) = 0;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP27
fun = @(x) 0.01*(x(1)-1)^2 + (x(2)-x(1)^2)^2;
nlcon = @(x) x(1)+x(3)^2 + 1;
cl = 0;
cu = 0;
x0 = [2;2;2];
fmin(i) = 0.04;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP28
fun = @(x) (x(1)+x(2))^2 + (x(2)+x(3))^2;
nlcon = @(x) x(1)+2*x(2)+3*x(3) - 1; %note linear
cl = 0;
cu = 0;
x0 = [-4;1;1];
fmin(i) = 0;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP29
fun = @(x) -x(1)*x(2)*x(3);
nlcon = @(x) -x(1)^2 - 2*x(2)^2 - 4*x(3)^2 + 48;
cl = 0;
cu = inf;
x0 = [1;1;1];
fmin(i) = -16*sqrt(2);

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP30
fun = @(x) x(1)^2 + x(2)^2 + x(3)^2;
nlcon = @(x) x(1)^2 + x(2)^2 - 1;
cl = 0;
cu = inf;
lb = [1;-10;-10];
ub = [10;10;10];
x0 = [1;1;1];
fmin(i) = 1;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,ub,nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP31
fun = @(x) 9*x(1)^2 + x(2)^2 + 9*x(3)^2;
nlcon = @(x) x(1)*x(2) - 1;
cl = 0;
cu = inf;
lb = [-10;1;-10];
ub = [10;10;1];
x0 = [1;1;1];
fmin(i) = 6;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,ub,nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP32
fun = @(x) (x(1) + 3*x(2) + x(3))^2 + 4*(x(1) - x(2))^2;
nlcon = @(x) [6*x(2) + 4*x(3) - x(1)^3 - 3;
              1 - x(1) - x(2) - x(3)];
cl = [0;0];
cu = [inf;0];
lb = [0;0;0];
x0 = [0.1;0.7;0.2];
fmin(i) = 1;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP33
fun = @(x) (x(1)-1)*(x(1)-2)*(x(1)-3)+x(3);
nlcon = @(x) [x(3)^2 - x(2)^2 - x(1)^2;
              x(1)^2 + x(2)^2 + x(3)^2 - 4];
cl = [0;0];
cu = [inf;inf];
lb = [0;0;0];
ub = [inf;inf;5];
x0 = [0;0;3];
fmin(i) = sqrt(2)-6;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,ub,nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP34
fun = @(x) -x(1);
nlcon = @(x) [x(2) - exp(x(1));
              x(3) - exp(x(2))];
cl = [0;0];
cu = [inf;inf];
lb = [0;0;0];
ub = [100;100;10];
x0 = [0;1.05;2.9];
fmin(i) = -log(log(10));

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,ub,nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP35
fun = @(x) 9 - 8*x(1) - 6*x(2) - 4*x(3) + 2*x(1)^2 + 2*x(2)^2 + x(3)^2 + 2*x(1)*x(2) + 2*x(1)*x(3);
nlcon = @(x) 3 - x(1) - x(2) - 2*x(3); 
cl = 0;
cu = inf;
lb = [0;0;0];
x0 = [0.5;0.5;0.5];
fmin(i) = 1/9; 

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP36
fun = @(x) -x(1)*x(2)*x(3);
nlcon = @(x) 72 - x(1) - 2*x(2) - 2*x(3);
cl = 0;
cu = inf;
lb = [0;0;0];
ub = [20;11;42];
x0 = [10;10;10];
fmin(i) = -3300;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,ub,nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP37
fun = @(x) -x(1)*x(2)*x(3);
nlcon = @(x) [72 - x(1) - 2*x(2) - 2*x(3);
              x(1) + 2*x(2) + 2*x(3)]; 
cl = [0;0];
cu = [inf;inf];
lb = [0;0;0];
ub = [42;42;42];
x0 = [10;10;10];
fmin(i) = -3456;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,ub,nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP38
fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2 + 90*(x(4)-x(3)^2)^2 + (1-x(3))^2 + 10.1*((x(2)-1)^2 + (x(4)-1)^2) + 19.8*(x(2)-1)*(x(4)-1);
lb = [-10;-10;-10;-10];
ub = [10;10;10;10];
x0 = [-3;-1;-3;-1];
fmin(i) = 0;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,ub,[],[],[],[],x0,opts);
i = i + 1;

%% NLP39
fun = @(x) -x(1);
nlcon = @(x) [x(2) - x(1)^3 - x(3)^2;
              x(1)^2 - x(2) - x(4)^2];
cl = [0;0];
cu = [0;0];
x0 = [2;2;2;2];
fmin(i) = -1;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP40
fun = @(x) -x(1)*x(2)*x(3)*x(4);
nlcon = @(x) [x(1)^3 + x(2)^2 - 1;
              x(1)^2*x(4) - x(3);
              x(4)^2 - x(2)];
cl = [0;0;0];
cu = [0;0;0];
x0 = [0.8;0.8;0.8;0.8];
fmin(i) = -0.25;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP41
fun = @(x) 2 - x(1)*x(2)*x(3);
nlcon = @(x) x(1) + 2*x(2) + 2*x(3) - x(4); 
cl = 0;
cu = 0;
lb = [0;0;0;0];
ub = [1;1;1;2];
x0 = [2;2;2;2];
fmin(i) = 52/27;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,ub,nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP42
fun = @(x) (x(1)-1)^2 + (x(2)-2)^2 + (x(3)-3)^2 + (x(4)-4)^2;
nlcon = @(x) [x(1) - 2;
              x(3)^2 + x(4)^2 - 2];
cl = [0;0];
cu = [0;0];
x0 = [1;1;1;1];
fmin(i) = 28-10*sqrt(2);

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP43
fun = @(x) x(1)^2 + x(2)^2 + 2*x(3)^2 + x(4)^2 - 5*x(1) - 5*x(2) - 21*x(3) + 7*x(4);
nlcon = @(x) [8 - x(1)^2 - x(2)^2 - x(3)^2 - x(4)^2 - x(1) + x(2) - x(3) + x(4);
              10 - x(1)^2 - 2*x(2)^2 - x(3)^2 - 2*x(4)^2 + x(1) + x(4);
              5 - 2*x(1)^2 - x(2)^2 - x(3)^2 - 2*x(1) + x(2) + x(4)];
cl = [0;0;0];
cu = [inf;inf;inf];
x0 = [0;0;0;0];
fmin(i) = -44;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP44
fun = @(x) x(1) - x(2) - x(3) - x(1)*x(3) + x(1)*x(4) + x(2)*x(3) - x(2)*x(4);
nlcon = @(x) [8 - x(1) - 2*x(2); 
              12 - 4*x(1) - x(2);
              12 - 3*x(1) - 4*x(2);
              8 - 2*x(3) - x(4);
              8 - x(3) - 2*x(4);
              5 - x(3) - x(4)];
cl = [0;0;0;0;0;0];
cu = [inf;inf;inf;inf;inf;inf];
lb = [0;0;0;0];
x0 = [0;0;0;0];
fmin(i) = -15;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP45
fun = @(x) 2 - 1/120*x(1)*x(2)*x(3)*x(4)*x(5);
lb = [0;0;0;0;0];
ub = [1;2;3;4;5];
x0 = [2;2;2;2;2];
fmin(i) = 1;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,ub,[],[],[],[],x0,opts);
i = i + 1;

%% NLP46 (won't run as trig)
% fun = @(x) (x(1)-x(2))^2+(x(3)-1)^2+(x(4)-1)^4+(x(5)-1)^6;
% nlcon = @(x) [x(1)^2*x(4) + sin(x(4)-x(5)) - 1;
%               x(2) + x(3)^4*x(4)^2 - 2];
% cl = [0;0];
% cu = [0;0];
% x0 = [0.5*sqrt(2);1.75;.5;2;2];
% fmin(i) = 0;
% 
% [~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);
% i = i + 1;

%% NLP47
fun = @(x) (x(1)-x(2))^2 + (x(2)-x(3))^3 + (x(3)-x(4))^4 + (x(4)-x(5))^4;
nlcon = @(x) [x(1) + x(2)^2 + x(3)^3 - 3;
              x(2) - x(3)^2 + x(4) - 1;
              x(1)*x(5) - 1];
cl = [0;0;0];
cu = [0;0;0];
x0 = [2;sqrt(2);-1;2*-sqrt(2);0.5];
fmin(i) = -0.026715;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP48
fun = @(x) (x(1)-1)^2 + (x(2)-x(3))^2 + (x(4)-x(5))^2;
nlcon = @(x) [x(1) + x(2) + x(3) + x(4) + x(5) - 5;
              x(3) - 2*(x(4) + x(5)) + 3];
cl = [0;0];
cu = [0;0];
x0 = [3;5;-3;2;-2];
fmin(i) = 0;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP49
clc
fun = @(x) (x(1)-x(2))^2 + (x(3)-1)^2 + (x(4)-1)^4 + (x(5)-1)^6;
nlcon = @(x) [x(1) + x(2) + x(3) + 4*x(4) - 7;
              x(3) + 5*x(5) - 6];
cl = [0;0];
cu = [0;0];
x0 = [10;7;2;-3;0.8];
fmin(i) = 0;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP50
fun = @(x) (x(1)-x(2))^2 + (x(2)-x(3))^2 + (x(3)-x(4))^4 + (x(4)-x(5))^2;
nlcon = @(x) [x(1) + 2*x(2) + 3*x(3) - 6;
              x(2) + 2*x(3) + 3*x(4) - 6;
              x(3) + 2*x(4) + 3*x(5) - 6];
cl = [0;0;0];
cu = [0;0;0];
x0 = [35;-31;11;5;-5];
fmin(i) = 0;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);
i = i + 1;

%% NLP51 Rastrigin's Function (won't work - trig)
% fun = @(x) 20 + x(1)^2 + x(2)^2 - 10*(cos(2*pi*x(1)) + cos(2*pi*x(2)));
% lb = [-5;-5];
% ub = [5;5];
% fmin(i) = 0;
% 
% [~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,ub);
% i = i + 1;

%% NLP52 Sawtooth XY (won't work - trig)
% fun = @(x) ((sin(x(1)) - sin(2*x(1))/2 + sin(3*x(1))/3 - sin(4*x(1))/4 + 4) * (x(1)^2/(x(1) + 1))) * (2 + cos(x(2)) + cos(2*x(2) - 0.5)/2);
% x0 = [100;-50];
% fmin(i) = 0;
% 
% [~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],[],[],[],[],x0,opts);
% i = i + 1;

%% MINLP1
clc
fun = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
A = [-1 1;
      1 1];
rl = [-inf;5];
ru = [-1;5];
lb = [0;0]; 
ub = [4;4];
xtype = 'IC';
x0 = [2;2];
fmin(i) = 4904;

[~,fval(i),ef(i)] = opti_scipnl(fun,A,rl,ru,lb,ub,[],[],[],xtype,x0);
i = i + 1;

%% MINLP2
fun = @(x) (x(1) - 5)^2 + x(2)^2 - 25;
nlcon = @(x) -x(1)^2 + x(2)-0.5;
cl = 0;
cu = inf;      
xtype = 'II';
x0 = [4.9;0.1];
fmin(i) = -5;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,xtype,x0);
i = i + 1;

%% MINLP3
fun = @(x) -x(1) - x(2) - x(3);
nlcon = @(x) [ (x(2) - 1./2.)*(x(2) - 1./2.) + (x(3) - 1./2.)*(x(3) - 1./2.);
                x(1) - x(2);
                x(1) + x(3) + x(4)];
cl = [-inf;-inf;-inf];            
cu = [1/4;0;2];
ub = [1;Inf;Inf;5];
lb = [0;0;0;0];
xtype = 'BCCI';
x0 = [0;0;0;0];
fmin(i) = -2.50001414;

[~,fval(i),ef(i)] = opti_scipnl(fun,[],[],[],lb,ub,nlcon,cl,cu,xtype,x0);
i = i + 1;

%% Final Accuracy Details
clc
len = length(fmin);
fprintf('\nSCIP - Interface [MI]NLP Testing\n\n');
fprintf(' No     SCIP Fmin      Saved Fmin    Rel Error   Exitflag  Check\n')
relerr = zeros(len,1);
for i = 1:len
    if(fmin(i))
        relerr(i) = abs((fmin(i)-fval(i))/fmin(i));
    else
        relerr(i) = abs(fmin(i)-fval(i));
    end
    if(relerr(i) > 1e-3), str = 'Check'; else str = []; end
    fprintf('%3d: %12.6f  | %12.6f  | %10.2g  |   %d    |%6s |\n',i,fval(i),fmin(i),relerr(i),ef(i),str);
end
