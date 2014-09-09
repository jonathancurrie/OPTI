function [prob,sol,fmin] = nlp_HS(varargin)
%NLP_HS  Return an OPTI NLP from the Hock-Schittkowski Collection
%
%   prob = nlp_HS(no) return a pre-built optiprob of an HS problem no.
%
%   [prob,sol,fmin] = nlp_HS(no) returns the optimum solution and function
%   eval at the optimum
%
%   no = nlp_HS() returns the number of problems available for testing

%   (C) 2011 Jonathan Currie (I2C2)

% Problem Description Key:
%
% e.g. OCD-Kr-s
%
%   O:  Objective Function
%       C = Constant
%       L = Linear
%       Q = Quadratic
%       S = Sum of Squares
%       P = Generalized Polynomial
%       G = General
%
%   C:  Constraints
%       U = Unconstrained
%       B = Bounds
%       L = Linear
%       Q = Quadratic
%       P = Generalized Polynomial
%       G = General
%
%   D:  Regularity (regular if 1st and 2nd derivatives exist)
%       R = Regular
%       I = Irregular
%
%   K:  Solution
%       T = Known
%       P = Not Known
%
%   r:  Order of Partial Derivatives
%       0 = derivatives not implemented
%       1 = first derivatives implemented
%
%   s:  Serial number within class (i.e. problem number)

 %#ok<*NBRAK>

%Check if just returning no problems
if(nargin < 1)
    prob = 45; sol = []; fmin = [];
    return;
else
    no = varargin{1};
end

%Big switch yard
switch(no)
    case 1 %PBR-T1-1 Betts[8]
        fun = @(x) 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
        grad = @(x)[[2*x(1)-400*x(1)*(-x(1)^2+x(2))-2];[-200*x(1)^2+200*x(2)]];
        lb = [-inf; -1.5];
        x0 = [-2; 1];        
        prob = optiprob('obj',fun,'grad',grad,'lb',lb,'x0',x0);            
        sol = [1;1];
        fmin = 0;
        
    case 2 %PBR-T1-2 Betts[8]
        fun = @(x) 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
        grad = @(x)[[2*x(1)-400*x(1)*(-x(1)^2+x(2))-2];[-200*x(1)^2+200*x(2)]];
        lb = [-inf; 1.5];
        x0 = [-2; 1];        
        prob = optiprob('obj',fun,'grad',grad,'lb',lb,'x0',x0);            
        sol = [1.224370749;1.5];
        fmin = 0.0504261879;
        
    case 3 %QBR-T1-1 Schuldt[56]
        fun = @(x) x(2) + 1e-5*(x(2)-x(1))^2;
        grad = @(x)[[0.00002*x(1)-0.00002*x(2)];[0.00002*x(2)-0.00002*x(1)+1]];
        lb = [-inf;0];
        x0 = [10;1];
        prob = optiprob('obj',fun,'grad',grad,'lb',lb,'x0',x0);
        sol = [0;0];
        fmin = 0;
        
    case 4 %PBR-T1-3 Asaadi[1]
        fun = @(x) 1/3*(x(1) + 1)^3 + x(2);
        grad = @(x)[[(x(1)+1)^2];[1]];
        lb = [1;0];
        x0 = [1.125;0.125];
        prob = optiprob('obj',fun,'grad',grad,'lb',lb,'x0',x0);
        sol = [1;0];
        fmin = 8/3;
        
    case 5 %GBR-T1-1 McCormick[41]
        fun = @(x) sin(x(1) + x(2)) + (x(1) - x(2))^2 - 1.5*x(1) + 2.5*x(2) + 1;
        grad = @(x)[[2*x(1)-2*x(2)+cos(x(1)+x(2))-1.5];[2*x(2)-2*x(1)+cos(x(1)+x(2))+2.5]];     
        lb = [-1.5;-3];
        ub = [4;3];
        x0 = [0;0];
        prob = optiprob('obj',fun,'grad',grad,'bounds',lb,ub,'x0',x0);
        sol = [-pi/3+0.5;-pi/3-0.5];
        fmin = -0.5*sqrt(3)-pi/3;
        
    case 6 %QQR-T1-1 Betts[8]
        fun = @(x) (1-x(1))^2;
        grad = @(x)[[2*x(1)-2];[0]]; 
        nlcon = @(x) 10*(x(2) - x(1)^2);
        nljac = @(x)[[-20*x(1),10]];
        nlrhs = 0;
        nle = 0;
        x0 = [-1.2;1];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'x0',x0);
        sol = [1;1];
        fmin = 0;
        
    case 7 %GPR-T1-1 Miele e.al.[44,45]
        fun = @(x) log(1+x(1)^2) - x(2);
        grad = @(x)[[(2*x(1))/(x(1)^2+1)];[-1]];
        nlcon = @(x) (1 + x(1)^2)^2 + x(2)^2 - 4;
        nljac = @(x)[[4*x(1)*(x(1)^2+1),2*x(2)]];
        nlrhs = 0;
        nle = 0;
        x0 = [2;2];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'x0',x0);
        sol = [0;sqrt(3)];
        fmin = -sqrt(3);
        
    case 8 %CQR-T1-1 Betts[8]
        fun = @(x) -1;
        grad = @(x)[[0];[0]];
        nlcon = @(x) [x(1)^2 + x(2)^2 - 25;
                      x(1)*x(2)-9];
        nljac = @(x)[[2*x(1),2*x(2)];[x(2),x(1)]];
        nlrhs = [0;0];
        nle = [0;0];
        x0 = [2;1];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'x0',x0);
        sol = [sqrt(0.5*(25+sqrt(301)));9/sqrt(0.5*(25+sqrt(301)))];
        fmin = -1;
        
    case 9 %GLR-T1-1 Miele e.al.[44]
        fun = @(x) sin(pi*x(1)/12)*cos(pi*x(2)/16);
        grad = @(x)[[(pi*cos((pi*x(1))/12)*cos((pi*x(2))/16))/12];[-(pi*sin((pi*x(1))/12)*sin((pi*x(2))/16))/16]];
        nlcon = @(x) 4*x(1)-3*x(2);
        nljac = @(x)[[4,-3]];
        nlrhs = 0;
        nle = 0;
        x0 = [0;0];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'x0',x0);
        sol = [-39;-52]; %multiple solutions...
        fmin = -.5;
        
    case 10 %LQR-T1-1 Biggs[10]
        fun = @(x) x(1) - x(2);
        grad = @(x)[[1];[-1]];
        nlcon = @(x) -3*x(1)^2 + 2*x(1)*x(2) - x(2)^2 + 1;
        nljac = @(x)[[2*x(2)-6*x(1),2*x(1)-2*x(2)]];
        nlrhs = 0;
        nle = 1;
        x0 = [-10;10];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'x0',x0);
        sol = [0,1]; 
        fmin = -1;
        
    case 11 %QQR-T1-2 Biggs[10]
        fun = @(x) (x(1) - 5)^2 + x(2)^2 - 25;
        grad = @(x)[[2*x(1)-10];[2*x(2)]];
        nlcon = @(x) -x(1)^2 + x(2);
        nljac = @(x)[[-2*x(1),1]];
        nlrhs = 0;
        nle = 1;
        x0 = [4.9;0.1];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'x0',x0);
        a = 7.5*sqrt(6) + sqrt(338.5);
        sol = [(a-1/a)/sqrt(6),(a^2-2+a^(-2))/6]; %seems wrong??
        fmin = -8.498464223;
        
    case 12 %QQR-T1-3 Mine e.al.[46]
        fun = @(x) 0.5*x(1)^2 + x(2)^2 - x(1)*x(2) - 7*x(1) - 7*x(2);
        grad = @(x)[[1.0*x(1)-x(2)-7];[2*x(2)-x(1)-7.0]];
        nlcon = @(x) 25-4*x(1)^2-x(2)^2;
        nljac = @(x)[[-8*x(1),-2*x(2)]];
        nlrhs = 0;
        nle = 1;
        x0 = [0;0];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'x0',x0);
        sol = [2;3];
        fmin = -30;
        
    case 13 %QPR-T1-1 Betts[8], Kuhn, Tucker[38]
        fun = @(x) (x(1)-2)^2 + x(2)^2;
        grad = @(x)[[2*x(1)-4];[2*x(2)]];
        nlcon = @(x)(1-x(1))^3 - x(2);
        nljac = @(x)[[-3*(x(1)-1)^2,-1]];
        nlrhs = 0;
        nle = 1;
        lb = [0;0];
        x0 = [-2;-2];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'lb',lb,'x0',x0);
        sol = [1;0];
        fmin = 0.994578532840774;
        
    case 14 %QQR-T1-4 Bracken, McCormick[13], Himmelblau[29]
        fun = @(x) (x(1)-2)^2 + (x(2)-1)^2;
        grad = @(x)[[2*x(1)-4];[2*x(2)-2]];
        nlcon = @(x) [-0.25*x(1)^2 - x(2)^2 + 1;
                      x(1) - 2*x(2) + 1];
        nljac = @(x)[[-0.5*x(1),-2*x(2)];[1,-2]];
        nlrhs = [0;0];
        nle = [1;0];
        x0 = [2;2];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'x0',x0);
        sol = [0.5*(sqrt(7)-1);0.25*(sqrt(7)+1)];
        fmin = 9-2.875*sqrt(7);
        
    case 15 %PQR-T1-1 Betts[8]
        fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
        grad = @(x)[[2*x(1)-400*x(1)*(-x(1)^2+x(2))-2];[-200*x(1)^2+200*x(2)]];
        nlcon = @(x) [x(1)*x(2)-1
                      x(1) + x(2)^2];
        nljac = @(x)[[x(2),x(1)];[1,2*x(2)]];
        nlrhs = [0;0];
        nle = [1;1];
        ub = [0.5;inf];
        x0 = [-2,1];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'ub',ub,'x0',x0);
        sol = [0.5;2];
        fmin = 306.5;
        
    case 16 %PQR-T1-2 Betts[8]
        fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
        grad = @(x)[[2*x(1)-400*x(1)*(-x(1)^2+x(2))-2,-200*x(1)^2+200*x(2)]];
        nlcon = @(x) [x(1) + x(2)^2
                      x(1)^2 + x(2)];
        nljac = @(x)[[1,2*x(2)];[2*x(1),1]];
        nlrhs = [0;0];
        nle = [1;1];
        lb = [-0.5;-inf];
        ub = [0.5;1];
        x0 = [-2,1];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'x0',x0);
        sol = [0.5;0.25];
        fmin = 0.25;
        
    case 17 %PQR-T1-3 Betts[8]
        fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
        grad = @(x)[[2*x(1)-400*x(1)*(-x(1)^2+x(2))-2];[-200*x(1)^2+200*x(2)]];
        nlcon = @(x) [x(2)^2 - x(1)
                      x(1)^2 - x(2)];
        nljac = @(x)[[-1,2*x(2)];[2*x(1),-1]];
        nlrhs = [0;0];
        nle = [1;1];
        lb = [-0.5;-inf];
        ub = [0.5;1];
        x0 = [-2,1];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'x0',x0);
        sol = [0;0];
        fmin = 1;
    
    case 18 %QQR-T1-5 Betts[8]
        fun = @(x) 0.01*x(1)^2 + x(2)^2;
        grad = @(x)[[0.02*x(1)];[2*x(2)]];
        nlcon = @(x) [x(1)*x(2) - 25;
                      x(1)^2 + x(2)^2 - 25];
        nljac = @(x)[[x(2),x(1)];[2*x(1),2*x(2)]];
        nlrhs = [0;0];
        nle = [1;1];
        lb = [2;0];
        ub = [50;50];
        x0 = [2;2];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'x0',x0);
        sol = [sqrt(250),sqrt(2.5)];
        fmin = 5;
        
    case 19 %PQR-T1-4 Betts[8], Gould[27]
        fun = @(x) (x(1)-10)^3 + (x(2)-20)^3;
        grad = @(x)[[3*(x(1)-10)^2];[3*(x(2)-20)^2]];
        nlcon = @(x) [(x(1)-5)^2 + (x(2)-5)^2 - 100;
                      -(x(2)-5)^2 - (x(1)-6)^2 + 82.81];
        nljac = @(x)[[2*x(1)-10,2*x(2)-10];[12.0-2*x(1),10.0-2*x(2)]];
        nlrhs = [0;0];
        nle = [1;1];
        lb = [13;0];
        ub = [100;100];
        x0 = [20.1;5.84];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'x0',x0);
        sol = [14.095;0.84296079];
        fmin = -6961.81381;
        
    case 20 %PQR-T1-5 Betts[8]
        fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
        grad = @(x)[[2*x(1)-400*x(1)*(-x(1)^2+x(2))-2];[-200*x(1)^2+200*x(2)]];
        nlcon = @(x) [x(1) + x(2)^2;
                      x(1)^2 + x(2);
                      x(1)^2 + x(2)^2 - 1];
        nljac = @(x)[[1,2*x(2)];[2*x(1),1];[2*x(1),2*x(2)]];
        nlrhs = [0;0;0];
        nle = [1;1;1];
        lb = [-0.5;-inf];
        ub = [0.5;inf];
        x0 = [-2;1];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'x0',x0);
        sol = [0.5;0.5*sqrt(3)];
        fmin = 81.5-25*sqrt(3);
        
    case 21 %QLR-T1-1 Betts[8]
        fun = @(x) 0.01*x(1)^2 + x(2)^2 - 100;
        grad = @(x)[[0.02*x(1)];[2*x(2)]];
        nlcon = @(x) 10*x(1) - x(2) - 10; %note linear
        nljac = @(x)[[10,-1]];
        nlrhs = 0;
        nle = 1;
        lb = [2;-50];
        ub = [50;50];
        x0 = [-1;-1];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'x0',x0);
        sol = [2;0];
        fmin = -99.96;
        
    case 22 %QQR-T1-6 Bracken, McCormick[13], Himmelblau[29], Sheela[57]
        fun = @(x) (x(1)-2)^2 + (x(2)-1)^2;
        grad = @(x)[[2*x(1)-4];[2*x(2)-2]];
        nlcon = @(x) [-x(1) - x(2) + 2; %note linear
                      -x(1)^2 + x(2)];
        nljac = @(x)[[-1,-1];[-2*x(1),1]];
        nlrhs = [0;0];
        nle = [1;1];
        x0 = [2;2];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'x0',x0);
        sol = [1;1];
        fmin = 1;
        
    case 23 %QQR-T1-7 Betts[8]
        fun = @(x) x(1)^2 + x(2)^2;
        grad = @(x)[[2*x(1)];[2*x(2)]];
        nlcon = @(x) [x(1) + x(2) - 1;      %note linear
                      x(1)^2 + x(2)^2 - 1;
                      9*x(1)^2 + x(2)^2 - 9;
                      x(1)^2 - x(2);
                      x(2)^2 - x(1)];
        nljac = @(x)[[1,1];[2*x(1),2*x(2)];[18*x(1),2*x(2)];[2*x(1),-1];[-1,2*x(2)]];
        nlrhs = [0;0;0;0;0];
        nle = [1;1;1;1;1];
        lb = [-50;-50];
        ub = [50;50];
        x0 = [3;1];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'x0',x0);
        sol = [1;1];
        fmin = 2;
        
    case 24 %PLR-T1-1 Betts[8], Box[12]
        fun = @(x) 1/(27*sqrt(3)) * ((x(1) - 3)^2 - 9)*x(2)^3;
        grad = @(x)[[(3^(1/2)*x(2)^3*(2*x(1)-6))/81];[(3^(1/2)*x(2)^2*((x(1)-3)^2-9))/27]];
        nlcon = @(x) [x(1)/sqrt(3) - x(2);
                      x(1) + sqrt(3)*x(2);
                      -x(1) - sqrt(3)*x(2) + 6];
        nljac = @(x)[[3^(1/2)/3,-1];[1,3^(1/2)];[-1,-3^(1/2)]];
        nlrhs = [0;0;0];
        nle = [1;1;1];
        lb = [0;0];
        x0 = [1;0.5];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'lb',lb,'x0',x0);
        sol = [3;sqrt(3)];
        fmin = -1;
        
    case 25 %SBR-T1-1 Holzman[32], Himmelblau[29]
        fun = @prob25;
        lb = [0.1;0;0];
        ub = [100;25.6;5];
        x0 = [100;12.5;3];
        prob = optiprob('obj',fun,'bounds',lb,ub,'x0',x0);
        sol = [50;25;1.5];
        fmin = 0;
        
    case 26 %PPR-T1-1 Huang, Aggerwal[34], Miele e.al.[43]
        fun = @(x) (x(1)-x(2))^2 + (x(2)-x(3))^4;
        grad = @(x)[[2*x(1)-2*x(2)];[2*x(2)-2*x(1)+4*(x(2)-x(3))^3];[-4*(x(2)-x(3))^3]];
        nlcon = @(x) (1+x(2)^2)*x(1) + x(3)^4 - 3;
        nljac = @(x)[[x(2)^2+1,2*x(1)*x(2),4*x(3)^3]];
        nlrhs = 0;
        nle = 0;
        x0 = [-2.6;2;2];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'x0',x0);
        sol = [1;1;1];
        fmin = 0;
        
    case 27 %PQR-T1-6 Miele e.al.[44,45]
        fun = @(x) 0.01*(x(1)-1)^2 + (x(2)-x(1)^2)^2;
        grad = @(x)[[0.02*x(1)-4*x(1)*(-x(1)^2+x(2))-0.02];[-2*x(1)^2+2*x(2)];[0]];
        nlcon = @(x) x(1)+x(3)^2 + 1;
        nljac = @(x)[[1,0,2*x(3)]];
        nlrhs = 0;
        nle = 0;
        x0 = [2;2;2];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'x0',x0);
        sol = [-1;1;0];
        fmin = 0.04;
        
    case 28 %QLR-T1-2 Huang, Aggerwal[34]
        fun = @(x) (x(1)+x(2))^2 + (x(2)+x(3))^2;
        grad = @(x)[[2*x(1)+2*x(2)];[2*x(1)+4*x(2)+2*x(3)];[2*x(2)+2*x(3)]];
        nlcon = @(x) x(1)+2*x(2)+3*x(3) - 1; %note linear
        nljac = @(x)[[1,2,3]];
        nlrhs = 0;
        nle = 0;
        x0 = [-4;1;1];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'x0',x0);
        sol = [0.5;-0.5;0.5];
        fmin = 0;
        
    case 29 %PQR-T1-7 Biggs[10]
        fun = @(x) -x(1)*x(2)*x(3);
        grad = @(x)[[-x(2)*x(3)];[-x(1)*x(3)];[-x(1)*x(2)]];
        nlcon = @(x) -x(1)^2 - 2*x(2)^2 - 4*x(3)^2 + 48;
        nljac = @(x)[[-2*x(1),-4*x(2),-8*x(3)]];
        nlrhs = 0;
        nle = 1;
        x0 = [1;1;1];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'x0',x0);
        sol = [-4;-2*sqrt(2);2]; %multiple solutions
        fmin = -16*sqrt(2);
        
    case 30 %QQR-T1-8 Betts[8]
        fun = @(x) x(1)^2 + x(2)^2 + x(3)^2;
        grad = @(x)[[2*x(1)];[2*x(2)];[2*x(3)]];
        nlcon = @(x) x(1)^2 + x(2)^2 - 1;
        nljac = @(x)[[2*x(1),2*x(2),0]];
        nlrhs = 0;
        nle = 1;
        lb = [1;-10;-10];
        ub = [10;10;10];
        x0 = [1;1;1];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'x0',x0);
        sol = [1;0;0];
        fmin = 1;
        
    case 31 %QQR-T1-9 Betts[8]
        fun = @(x) 9*x(1)^2 + x(2)^2 + 9*x(3)^2;
        grad = @(x)[[18*x(1)];[2*x(2)];[18*x(3)]];
        nlcon = @(x) x(1)*x(2) - 1;
        nljac = @(x)[[x(2),x(1),0]];
        nlrhs = 0;
        nle = 1;
        lb = [-10;1;-10];
        ub = [10;10;1];
        x0 = [1;1;1];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'x0',x0);
        sol = [1/sqrt(3);sqrt(3);0];
        fmin = 6;
        
    case 32 %QPR-T1-2 Evtushenko[25]
        fun = @(x) (x(1) + 3*x(2) + x(3))^2 + 4*(x(1) - x(2))^2;
        grad = @(x)[[10*x(1)-2*x(2)+2*x(3)];[26*x(2)-2*x(1)+6*x(3)];[2*x(1)+6*x(2)+2*x(3)]];
        nlcon = @(x) [6*x(2) + 4*x(3) - x(1)^3 - 3;
                      1 - x(1) - x(2) - x(3)]; %note linear
        nljac = @(x)[[-3*x(1)^2,6,4];[-1,-1,-1]];
        nlrhs = [0;0];
        nle = [1;0];
        lb = [0;0;0];
        x0 = [0.1;0.7;0.2];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'lb',lb,'x0',x0);
        sol = [0;0;1];
        fmin = 1;
        
    case 33 %PQR-T1-8 Beltrami[6], Hartman[28]
        fun = @(x) (x(1)-1)*(x(1)-2)*(x(1)-3)+x(3);
        grad = @(x)[[(x(1)-1)*(x(1)-2)+(x(1)-1)*(x(1)-3)+(x(1)-2)*(x(1)-3)];[0];[1]];
        nlcon = @(x) [x(3)^2 - x(2)^2 - x(1)^2;
                      x(1)^2 + x(2)^2 + x(3)^2 - 4];
        nljac = @(x)[[-2*x(1),-2*x(2),2*x(3)];[2*x(1),2*x(2),2*x(3)]];
        nlrhs = [0;0];
        nle = [1;1];
        lb = [0;0;0];
        ub = [inf;inf;5];
        x0 = [0;0;3];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'x0',x0);
        sol = [0;sqrt(2);sqrt(2)];
        fmin = sqrt(2)-6;
        
    case 34 %LGR-T1-1 Eckhardt[24]
        fun = @(x) -x(1);
        grad = @(x)[[-1];[0];[0]];
        nlcon = @(x) [x(2) - exp(x(1));
                      x(3) - exp(x(2))];
        nljac = @(x)[[-exp(x(1)),1,0];[0,-exp(x(2)),1]];
        nlrhs = [0;0];
        nle = [1;1];
        lb = [0;0;0];
        ub = [100;100;10];
        x0 = [0;1.05;2.9];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'x0',x0);
        sol = [log(log(10));log(10);10];
        fmin = -log(log(10));
        
    case 35 %QLR-T1-3 Asaadi[1], Charalambous[18], Dimitru[23], Sheela[57]
        fun = @(x) 9 - 8*x(1) - 6*x(2) - 4*x(3) + 2*x(1)^2 + 2*x(2)^2 + x(3)^2 + 2*x(1)*x(2) + 2*x(1)*x(3);
        grad = @(x)[[4*x(1)+2*x(2)+2*x(3)-8];[2*x(1)+4*x(2)-6];[2*x(1)+2*x(3)-4]];
        nlcon = @(x) 3 - x(1) - x(2) - 2*x(3); %note linear
        nljac = @(x)[[-1,-1,-2]];
        nlrhs = 0;
        nle = 1;
        lb = [0;0;0];
        x0 = [0.5;0.5;0.5];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'lb',lb,'x0',x0);
        sol = [4/3;7/9;4/9];
        fmin = 1/9; 
        
    case 36 %PLR-T1-2 Biggs[10]
        fun = @(x) -x(1)*x(2)*x(3);
        grad = @(x)[[-x(2)*x(3)];[-x(1)*x(3)];[-x(1)*x(2)]];
        nlcon = @(x) 72 - x(1) - 2*x(2) - 2*x(3); %note linear
        nljac = @(x)[[-1,-2,-2]];
        nlrhs = 0;
        nle = 1;
        lb = [0;0;0];
        ub = [20;11;42];
        x0 = [10;10;10];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'x0',x0);
        sol = [20;11;15];
        fmin = -3300;
        
    case 37 %PLR-T1-3 Betts[8], Box[12]
        fun = @(x) -x(1)*x(2)*x(3);
        grad = @(x)[[-x(2)*x(3)];[-x(1)*x(3)];[-x(1)*x(2)]];
        nlcon = @(x) [72 - x(1) - 2*x(2) - 2*x(3); %note linear
                      x(1) + 2*x(2) + 2*x(3)];  %note linear
        nljac = @(x)[[-1,-2,-2];[1,2,2]];
        nlrhs = [0;0];
        nle = [1;1];
        lb = [0;0;0];
        ub = [42;42;42];
        x0 = [10;10;10];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'x0',x0);
        sol = [24;12;12];
        fmin = -3456;   
        
    case 38 %PBR-T1-4 Colville[20], Himmelblau[29]
        fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2 + 90*(x(4)-x(3)^2)^2 + (1-x(3))^2 + 10.1*((x(2)-1)^2 + (x(4)-1)^2) + 19.8*(x(2)-1)*(x(4)-1);
        grad = @(x)[[2*x(1)-400*x(1)*(-x(1)^2+x(2))-2.0];[-200*x(1)^2+220.2*x(2)+19.8*x(4)-40.0];[2*x(3)-360*x(3)*(-x(3)^2+x(4))-2.0];[-180*x(3)^2+19.8*x(2)+200.2*x(4)-40.0]];
        lb = [-10;-10;-10;-10];
        ub = [10;10;10;10];
        x0 = [-3;-1;-3;-1];
        prob = optiprob('obj',fun,'grad',grad,'bounds',lb,ub,'x0',x0);
        sol = [1;1;1;1];
        fmin = 0;
        
    case 39 %LPR-T1-1 Miele e.al.[44,45]
        fun = @(x) -x(1);
        grad = @(x)[[-1];[0];[0];[0]];
        nlcon = @(x) [x(2) - x(1)^3 - x(3)^2;
                      x(1)^2 - x(2) - x(4)^2];
        nljac = @(x)[[-3*x(1)^2,1,-2*x(3),0];[2*x(1),-1,0,-2*x(4)]];
        nlrhs = [0;0];
        nle = [0;0];
        x0 = [2;2;2;2];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'x0',x0);
        sol = [1;1;0;0];
        fmin = -1;
        
    case 40 %PPR-T1-2 Beltrami[6], Indusi[35]
        fun = @(x) -x(1)*x(2)*x(3)*x(4);
        grad = @(x)[[-x(2)*x(3)*x(4)];[-x(1)*x(3)*x(4)];[-x(1)*x(2)*x(4)];[-x(1)*x(2)*x(3)]];
        nlcon = @(x) [x(1)^3 + x(2)^2 - 1;
                      x(1)^2*x(4) - x(3);
                      x(4)^2 - x(2)];
        nljac = @(x)[[3*x(1)^2,2*x(2),0,0];[2*x(1)*x(4),0,-1,x(1)^2];[0,-1,0,2*x(4)]];
        nlrhs = [0;0;0];
        nle = [0;0;0];
        x0 = [0.8;0.8;0.8;0.8];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'x0',x0);
        sol = [2^(-1/3);2^(2*-1/4);(-1)^2*2^(-11/12);(-1)^2*2^(-1/4)]; %multiple solutions
        fmin = -0.25;
        
    case 41 %PLR-T1-4 Betts[8], Miele e.al.[42]
        fun = @(x) 2 - x(1)*x(2)*x(3);
        grad = @(x)[[-x(2)*x(3)];[-x(1)*x(3)];[-x(1)*x(2)];[0]];
        nlcon = @(x) x(1) + 2*x(2) + 2*x(3) - x(4); %note linear
        nljac = @(x)[[1,2,2,-1]];
        nlrhs = 0;
        nle = 0;
        lb = [0;0;0;0];
        ub = [1;1;1;2];
        x0 = [2;2;2;2];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'x0',x0);
        sol = [2/3;1/3;1/3;2];
        fmin = 52/27;
        
    case 42 %QQR-T1-10 Brusch[14]
        fun = @(x) (x(1)-1)^2 + (x(2)-2)^2 + (x(3)-3)^2 + (x(4)-4)^2;
        grad = @(x)[[2*x(1)-2];[2*x(2)-4];[2*x(3)-6];[2*x(4)-8]];
        nlcon = @(x) [x(1) - 2; %note linear
                      x(3)^2 + x(4)^2 - 2];
        nljac = @(x)[[1,0,0,0];[0,0,2*x(3),2*x(4)]];
        nlrhs = [0;0];
        nle = [0;0];
        x0 = [1;1;1;1];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'x0',x0);
        sol = [2;2;0.6*sqrt(2);0.8*sqrt(2)];
        fmin = 28-10*sqrt(2);
        
    case 43 %QQR-T1-11 Betts[18], Charalambous[18], Gould[27], Sheela[57]
        fun = @(x) x(1)^2 + x(2)^2 + 2*x(3)^2 + x(4)^2 - 5*x(1) - 5*x(2) - 21*x(3) + 7*x(4);
        grad = @(x)[[2*x(1)-5];[2*x(2)-5];[4*x(3)-21];[2*x(4)+7]];
        nlcon = @(x) [8 - x(1)^2 - x(2)^2 - x(3)^2 - x(4)^2 - x(1) + x(2) - x(3) + x(4);
                      10 - x(1)^2 - 2*x(2)^2 - x(3)^2 - 2*x(4)^2 + x(1) + x(4);
                      5 - 2*x(1)^2 - x(2)^2 - x(3)^2 - 2*x(1) + x(2) + x(4)];
        nljac = @(x)[[-2*x(1)-1,1-2*x(2),-2*x(3)-1,1-2*x(4)];[1-2*x(1),-4*x(2),-2*x(3),1-4*x(4)];[-4*x(1)-2,1-2*x(2),-2*x(3),1]];
        nlrhs = [0;0;0];
        nle = [1;1;1];
        x0 = [0;0;0;0];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'x0',x0);
        sol = [0;1;2;-1];
        fmin = -44;
        
    case 44 %QLR-T1-4 Konno[37]
        fun = @(x) x(1) - x(2) - x(3) - x(1)*x(3) + x(1)*x(4) + x(2)*x(3) - x(2)*x(4);
        grad = @(x)[[x(4)-x(3)+1];[x(3)-x(4)-1];[x(2)-x(1)-1];[x(1)-x(2)]];
        nlcon = @(x) [8 - x(1) - 2*x(2); %note ALL linear
                      12 - 4*x(1) - x(2);
                      12 - 3*x(1) - 4*x(2);
                      8 - 2*x(3) - x(4);
                      8 - x(3) - 2*x(4);
                      5 - x(3) - x(4)];
        nljac = @(x)[[-1,-2,0,0];[-4,-1,0,0];[-3,-4,0,0];[0,0,-2,-1];[0,0,-1,-2];[0,0,-1,-1]];
        nlrhs = [0;0;0;0;0;0];
        nle = [1;1;1;1;1;1];
        lb = [0;0;0;0];
        x0 = [0;0;0;0];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'lb',lb,'x0',x0);
        sol = [0;3;0;4];
        fmin = -15;
        
    case 45 %PBR-T1-5 Betts[8], Miele e.al.[42]
        fun = @(x) 2 - 1/120*x(1)*x(2)*x(3)*x(4)*x(5);
        grad = @(x)[[-(x(2)*x(3)*x(4)*x(5))/120];[-(x(1)*x(3)*x(4)*x(5))/120];[-(x(1)*x(2)*x(4)*x(5))/120];[-(x(1)*x(2)*x(3)*x(5))/120];[-(x(1)*x(2)*x(3)*x(4))/120]];
        lb = [0;0;0;0;0];
        ub = [1;2;3;4;5];
        x0 = [2;2;2;2;2];
        prob = optiprob('obj',fun,'grad',grad,'bounds',lb,ub,'x0',x0);
        sol = [1;2;3;4;5];
        fmin = 1;
        
    case 46
        fun = @(x) (x(1)-x(2))^2+(x(3)-1)^2+(x(4)-1)^4+(x(5)-1)^6;
        grad = @(x)[[2*x(1)-2*x(2)];[2*x(2)-2*x(1)];[2*x(3)-2];[4*(x(4)-1)^3];[6*(x(5)-1)^5]];
        nlcon = @(x) [x(1)^2*x(4) + sin(x(4)-x(5)) - 1;
                      x(2) + x(3)^4*x(4)^2 - 2];
        nljac = @(x)[[2*x(1)*x(4),0,0,x(1)^2+cos(x(4)-x(5)),-cos(x(4)-x(5))];[0,1,4*x(3)^3*x(4)^2,2*x(3)^4*x(4),0]];
        nlrhs = [0;0];
        nle = [0;0];
        x0 = [0.5*sqrt(2);1.75;.5;2;2];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'x0',x0);
        sol = [1;1;1;1;1];
        fmin = 0;
        
    case 47 %PPR-T1-3 Huang, Aggerwal[34], Miele e.al.[43]
        fun = @(x) (x(1)-x(2))^2 + (x(2)-x(3))^3 + (x(3)-x(4))^4 + (x(4)-x(5))^4;
        grad = @(x)[[2*x(1)-2*x(2)];[2*x(2)-2*x(1)+3*(x(2)-x(3))^2];[4*(x(3)-x(4))^3-3*(x(2)-x(3))^2];[4*(x(4)-x(5))^3-4*(x(3)-x(4))^3];[-4*(x(4)-x(5))^3]];
        nlcon = @(x) [x(1) + x(2)^2 + x(3)^3 - 3;
                      x(2) - x(3)^2 + x(4) - 1;
                      x(1)*x(5) - 1];
        nljac = @(x)[[1,2*x(2),3*x(3)^2,0,0];[0,1,-2*x(3),1,0];[x(5),0,0,0,x(1)]];
        nlrhs = [0;0;0];
        nle = [0;0;0];
        x0 = [2;sqrt(2);-1;2*-sqrt(2);0.5];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'x0',x0);
        sol = [1;1;1;1;1];
        fmin = 0;
        
    case 48 %QLR-T1-5 Huang, Aggerwal[34], Miele e.al.[43]
        fun = @(x) (x(1)-1)^2 + (x(2)-x(3))^2 + (x(4)-x(5))^2;
        grad = @(x)[[2*x(1)-2];[2*x(2)-2*x(3)];[2*x(3)-2*x(2)];[2*x(4)-2*x(5)];[2*x(5)-2*x(4)]];
        nlcon = @(x) [x(1) + x(2) + x(3) + x(4) + x(5) - 5;
                      x(3) - 2*(x(4) + x(5)) + 3];
        nljac = @(x)[[1,1,1,1,1];[0,0,1,-2,-2]];
        nlrhs = [0;0];
        nle = [0;0];
        x0 = [3;5;-3;2;-2];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'x0',x0);
        sol = [1;1;1;1;1];
        fmin = 0;
        
    case 49 %PLR-T1-5 Huang, Aggerwal[34]
        fun = @(x) (x(1)-x(2))^2 + (x(3)-1)^2 + (x(4)-1)^4 + (x(5)-1)^6;
        grad = @(x)[[2*x(1)-2*x(2)];[2*x(2)-2*x(1)];[2*x(3)-2];[4*(x(4)-1)^3];[6*(x(5)-1)^5]];
        nlcon = @(x) [x(1) + x(2) + x(3) + 4*x(4) - 7;
                      x(3) + 5*x(5) - 6];
        nljac = @(x)[[1,1,1,4,0];[0,0,1,0,5]];
        nlrhs = [0;0];
        nle = [0;0];
        x0 = [10;7;2;-3;0.8];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'x0',x0);
        sol = [1;1;1;1;1];
        fmin = 0;
        
    case 50 %PLR-T1-6 Huang, Aggerwal[34]
        fun = @(x) (x(1)-x(2))^2 + (x(2)-x(3))^2 + (x(3)-x(4))^4 + (x(4)-x(5))^2;
        grad = @(x)[[2*x(1)-2*x(2)];[4*x(2)-2*x(1)-2*x(3)];[2*x(3)-2*x(2)+4*(x(3)-x(4))^3];[2*x(4)-2*x(5)-4*(x(3)-x(4))^3];[2*x(5)-2*x(4)]];
        nlcon = @(x) [x(1) + 2*x(2) + 3*x(3) - 6;
                      x(2) + 2*x(3) + 3*x(4) - 6;
                      x(3) + 2*x(4) + 3*x(5) - 6];
        nljac = @(x)[[1,2,3,0,0];[0,1,2,3,0];[0,0,1,2,3]];
        nlrhs = [0;0;0];
        nle = [0;0;0];
        x0 = [35;-31;11;5;-5];
        prob = optiprob('obj',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'x0',x0);
        sol = [1;1;1;1;1];
        fmin = 0;
        
        
    otherwise
        error('Unknown or unimplemented problem');
end




%Functions not easily described as anonymous functions
function f = prob25(x)
f = 0;
for i = 1:99
    u = 25 + (-50*log(0.01*i))^(2/3);
    fi = -0.01*i + exp(-1/x(1) * (u - x(2))^x(3));
    f = f + fi^2;
end




        