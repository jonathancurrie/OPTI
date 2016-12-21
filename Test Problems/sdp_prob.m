function [prob,sol,fmin] = sdp_prob(varargin)
%SDP_PROB  Return an OPTI SDP 
%
%   prob = sdp_prob(no) return a pre-built optiprob of a saved SDP.
%
%   [prob,sol,fmin] = sdp_prob(no) returns the optimum solution and 
%   function eval at the optimum
%
%   no = sdp_prob() returns the number of problems available for testing.

%   (C) 2011 Jonathan Currie (IPL)

%Check if just returning no problems
if(nargin < 1)
    prob = 5; sol = []; fmin = [];
    return;
else
    no = varargin{1};
end          

%Big switch yard
switch(no)
    case 1 
        f = 1;
        A = eye(2);
        C = -[0 sqrt(2); sqrt(2) 0];
        sdp = sparse([C(:) A(:)]);
        prob = optiprob('f',f,'sdcone',sdp);            
        sol = sqrt(2);
        fmin = sqrt(2);
        
    case 2 
        f = [1 1]';
        C = -[0 2; 2 0];
        A0 = [1 0; 0 0];
        A1 = [0 0; 0 1];
        sdp = sparse([C(:) A0(:) A1(:)]);
        lb = [0;0];
        ub = [10;10];
        prob = optiprob('f',f,'sdcone',sdp,'bounds',lb,ub);             
        sol = [2;2];
        fmin = 4;
        
    case 3 
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
        prob = optiprob('f',f,'sdcone',sdp);          
        sol = [1.2;1.1;0.1;1.1];
        fmin = 1.2;       
        
    case 4 
        f = [1 1 1]';
        C = -[0 1 2; 1 0 3; 2 3 100];
        A0 = [1 0 0; 0 0 0; 0 0 0];
        A1 = [0 0 0; 0 1 0; 0 0 0];
        A2 = zeros(3);
        sdp = sparse([C(:) A0(:) A1(:) A2(:)]);
        lb = [10;0;0];
        ub = [1000;1000;1000];
        prob = optiprob('f',f,'sdcone',sdp,'bounds',lb,ub);             
        sol = [9.99999907896292;0.178713937638058;-9.21054037686845e-07];
        fmin = 10.1787148490962;
        
    case 5
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
        prob = optiprob('f',f,'sdcone',sdp);          
        sol = [-0.749999999967481;-0.999999999573634];
        fmin = 2.74999999779937;
        
    otherwise
        error('Problem not available or not implemented yet');
end