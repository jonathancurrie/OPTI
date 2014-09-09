%% Test Problems for SEDUMI Interface
clc
clear

%Check SEDUMI is available
if(~checkSolver('sedumi',0))
    fprintf(2,'\n\nSEDUMI IS NOT AVAILABLE, NO SEDUMI TESTS CAN BE RUN\n\n');
    return;
end

files = {'sdp_truss1.mat','sdp_truss2.mat','sdp_truss3.mat','sdp_truss4.mat','sdp_hinf01.mat',...
         'sdp_hinf02.mat','sdp_hinf03.mat','sdp_hinf04.mat','sdp_arch0.mat','sdp_arch2.mat',...
         'sdp_arch4.mat','sdp_qap05.mat','sdp_qap06.mat','sdp_qap07.mat','sdp_qap08.mat',...
         'sdp_mcp100.mat','sdp_mcp124-1.mat','sdp_mcp124-2.mat','sdp_mcp124-3.mat','sdp_mcp124-4.mat',...
         'sdp_mcp250-1.mat','sdp_mcp250-2.mat','sdp_mcp250-3.mat','sdp_mcp250-4.mat'};
     
len = length(files);
cfval = zeros(len,1);
sfval = zeros(len,1);
optsC = optiset('solver','csdp','maxtime',15);
optsS = optiset('solver','sedumi','maxtime',15);

fprintf('SEDUMI SDP Checking\n');
%For each problem
for i = 1:len
    %Read In Problem
    prob = sdpRead(files{i});
    %Solve
    fprintf('Running Problem %s....',files{i});
    [~,cfval(i)] = solve(opti(prob,optsC));
    [~,sfval(i)] = solve(opti(prob,optsS));
    fprintf(' Done\n');
end

%Display fvals
display([cfval sfval])
if(norm(cfval-sfval) > 1e-1) 
    error('Failed SeDuMi Checking');
else
    fprintf('Passed SeDuMi Checking\n');
end

% %% Complex Checking
% % Create a positive (semi)definite matrix:
% R=randn(5)+j*randn(5); R=R*R'; 
% % Make sure it's Hermitian:
% R=(R+R')/2;
% 
% % Primal formulation
% clear At b c K
% c=vec(eye(5));
% At=vec(R);
% b=1;
% K.s=5;
% K.scomplex=1;
% 
% [x,y,info]=sedumi(At,b,c,K);
% 
% % Normed eigenvector:
% e=x(1:5)/norm(x(1:5));
% % Verify the eigenvalue
% [real(e'*R*e), max(eig(R))]

