%% TESTING SDPA READING AND WRITING ROUTINES
clear all
clc
errtol = 1e-1; %LOW tol as non-deterministic with CSDP OpenMP

%% SDPA Files To Read, Solve, Write, Read, Solve (dense and sparse)
files = {'sdp1.dat-s','sdp2.dat-s','sdp3.dat-s','sdp4.dat-s','truss1.dat-s','truss2.dat-s','truss3.dat-s','truss4.dat-s',...
         'sdpademo.dat','sdpademo.dat-s','sdpademo2.dat','hinf1.dat-s','hinf2.dat-s','hinf3.dat-s','hinf4.dat-s'};

len = length(files);
rfval = zeros(len,1);
wSfval = zeros(len,1);
wDfval = zeros(len,1);
opts = optiset('solver','csdp','maxtime',5);

fprintf('SDPA SDP Sparse+Dense Read+Write Checking\n');
%For each problem
for i = 1:len
    %Read In Problem
    prob = sdpRead(files{i});
    %Solve
    [~,rfval(i)] = solve(opti(prob,opts));
    %Write problem to new sparse and dense file
    name = regexp(files{i},'\.','split');
    fprintf('Checking: %s [%d Linear, %d SD]\n',name{1},size(prob.A,1),length(prob.sdcone));
    spname = [name{1} 'chk.dat-s'];
    dnname = [name{1} 'chk.dat'];
    sdpWrite(prob,spname);
    sdpWrite(prob,dnname);
    %Read the problem back from the written files
    probSP = sdpRead(spname);
    probDN = sdpRead(dnname);
    %Resolve each
    [~,wSfval(i)] = solve(opti(probSP,opts));
    [~,wDfval(i)] = solve(opti(probDN,opts));
    %Delete intermediate files
    delete(spname,dnname);
end

%Display fvals
display([rfval wSfval wDfval])
if(norm(rfval-wSfval) > errtol) 
    error('Failed SDPA Sparse Read/Write Check');
elseif(norm(rfval-wDfval) > errtol)
    error('Failed SDPA Dense Read/Write Check');
end

%% Check SDPA LP Problems (including ones with Infs)
clc
len = lp_prob;
rfval = zeros(len,1);
wSfval = zeros(len,1);
wDfval = zeros(len,1);
opts = optiset('solver','csdp','maxtime',5);

fprintf('SDPA LP Sparse+Dense Read+Write Checking\n');
%For each problem
for i = 1:len
    %Read In Problem
    prob = lp_prob(i);
    %Solve
    [~,rfval(i)] = solve(opti(prob,opts));
    %Write problem to new sparse and dense file
    name = prob.Name;
    fprintf('Checking: %s [%d Linear, %d SD]\n',name,size(prob.A,1),length(prob.sdcone));
    spname = [name 'chk.dat-s'];
    dnname = [name 'chk.dat'];
    sdpWrite(prob,spname);
    sdpWrite(prob,dnname);
    %Read the problem back from the written files
    probSP = sdpRead(spname);
    probDN = sdpRead(dnname);
    %Resolve each
    [~,wSfval(i)] = solve(opti(probSP,opts));
    [~,wDfval(i)] = solve(opti(probDN,opts));
    %Delete intermediate files
    delete(spname,dnname);
end

%Display fvals
display([rfval wSfval wDfval])
if(norm(rfval-wSfval) > errtol) %LOW tol until Brian corrects repeat problem error
    error('Failed SDPA (LP) Sparse Read/Write Check');
elseif(norm(rfval-wDfval) > errtol)
    error('Failed SDPA (LP) Dense Read/Write Check');
end


%% SDPA Files To Read, Solve, Write, Read, Solve (just sparse)
clc
files = {'arch0.dat-s','arch2.dat-s','arch4.dat-s','qap5.dat-s','qap6.dat-s','qap7.dat-s','qap8.dat-s',...
         'mcp100.dat-s','mcp124-1.dat-s','mcp124-2.dat-s','mcp124-3.dat-s','mcp124-4.dat-s',...
         'mcp250-1.dat-s','mcp250-2.dat-s','mcp250-3.dat-s','mcp250-4.dat-s'};

len = length(files);
rfval = zeros(len,1);
wSfval = zeros(len,1);
opts = optiset('solver','csdp','maxtime',15);

fprintf('SDPA Large SDP Sparse Read+Write Checking\n');
%For each problem
for i = 1:len
    %Read In Problem
    prob = sdpRead(files{i});
    %Solve
    [~,rfval(i)] = solve(opti(prob,opts));
    %Write problem to new sparse and dense file
    name = regexp(files{i},'\.','split');
    fprintf('Checking: %s [%d Linear, %d SD]\n',name{1},size(prob.A,1),length(prob.sdcone));
    spname = [name{1} 'chk.dat-s'];
    sdpWrite(prob,spname);
    %Read the problem back from the written files
    probSP = sdpRead(spname);
    %Resolve each
    [~,wSfval(i)] = solve(opti(probSP,opts));
    %Delete intermediate files
    delete(spname);
end

%Display fvals
display([rfval wSfval])
if(norm(rfval-wSfval) > errtol) %LOW tol until Brian corrects repeat problem error
    error('Failed SDPA Sparse Large Read/Write Check');
end


%% SEDUMI TESTING
clc
files = {'sdp_truss1.mat','sdp_truss2.mat','sdp_truss3.mat','sdp_truss4.mat','sdp_hinf01.mat',...
         'sdp_hinf02.mat','sdp_hinf03.mat','sdp_hinf04.mat','sdp_arch0.mat','sdp_arch2.mat',...
         'sdp_arch4.mat','sdp_qap05.mat','sdp_qap06.mat','sdp_qap07.mat','sdp_qap08.mat',...
         'sdp_mcp100.mat','sdp_mcp124-1.mat','sdp_mcp124-2.mat','sdp_mcp124-3.mat','sdp_mcp124-4.mat',...
         'sdp_mcp250-1.mat','sdp_mcp250-2.mat','sdp_mcp250-3.mat','sdp_mcp250-4.mat'};
     
len = length(files);
rfval = zeros(len,1);
wSfval = zeros(len,1);
opts = optiset('solver','csdp','maxtime',15);

fprintf('SEDUMI SDP Sparse Read+Write Checking\n');
%For each problem
for i = 1:len
    %Read In Problem
    prob = sdpRead(files{i});
    %Solve
    [~,rfval(i)] = solve(opti(prob,opts));
    %Write problem to new sparse and dense file
    name = regexp(files{i},'\.','split');
    fprintf('Checking: %s [%d Linear, %d SD]\n',name{1},prob.sdcone.K.l,length(prob.sdcone.K.s));
    spname = [name{1} 'chk.mat'];
    sdpWrite(prob,spname);
    %Read the problem back from the written files
    probSP = sdpRead(spname);
    %Resolve each
    [~,wSfval(i)] = solve(opti(probSP,opts));
    %Delete intermediate files
    delete(spname);
end

%Display fvals
display([rfval wSfval])
if(norm(rfval-wSfval) > errtol) %LOW tol until Brian corrects repeat problem error
    error('Failed SeDuMi Sparse Large Read/Write Check');
end

%% Check File Reading and Writing against each type (checks conversions)
clc
files = {'arch0','arch2','arch4','truss1','truss2','truss3','truss4','mcp100',...
         'mcp124-1','mcp124-2','mcp124-3','mcp124-4','mcp250-1','mcp250-2',...
         'mcp250-3','mcp250-4'};
     
len = length(files);
SDMIRval = zeros(len,1);
SDPARval = zeros(len,1);
SDMIWval = zeros(len,1);
SDPAWval = zeros(len,1);
opts = optiset('solver','csdp','maxtime',15);

fprintf('SEDUMI vs SDPA SDP Sparse Read+Write Checking\n');
%For each problem
for i = 1:len
    %Read In Problem
    prob1 = sdpRead([files{i} '.dat-s']);
    prob2 = sdpRead(['sdp_' files{i} '.mat']);
    %Solve
    fprintf('Running: %s [%d Linear, %d SD]\n',files{i},prob2.sdcone.K.l,length(prob2.sdcone.K.s));
    [~,SDPARval(i)] = solve(opti(prob1,opts));
    [~,SDMIRval(i)] = solve(opti(prob2,opts));
    %Write (swopping here)
    sdpWrite(prob1,[files{i} 'chk.mat']);
    sdpWrite(prob2,[files{i} 'chk.dat-s']);
    %Read and re-solve 
    prob1b = sdpRead([files{i} 'chk.dat-s']);
    prob2b = sdpRead([files{i} 'chk.mat']);
    [~,SDPAWval(i)] = solve(opti(prob1b,opts));
    [~,SDMIWval(i)] = solve(opti(prob2b,opts));
    %Delete temp files
    delete([files{i} 'chk.mat'],[files{i} 'chk.dat-s']);
end

%Display fvals
display([SDMIRval SDPARval SDMIWval SDPAWval])
if(norm(SDMIRval-SDPARval) > errtol || norm(SDMIWval-SDPAWval) > errtol || norm(SDMIRval-SDMIWval) > errtol) 
    error('Failed SDPA + SeDuMi Sparse Large Read + Write Check');
end     

%% Check SeDuMi LPS
clc
len = lp_prob;
rfval = zeros(len,1);
wSfval = zeros(len,1);
opts = optiset('solver','csdp','maxtime',5);

fprintf('SEDUMI LP Sparse Read+Write Checking\n');
%For each problem
for i = 1:len
    %Read In Problem
    prob = lp_prob(i);
    %Solve
    [~,rfval(i)] = solve(opti(prob,opts));
    %Write problem to new sparse and dense file
    name = prob.Name;
    fprintf('Checking: %s [%d Linear, %d SD]\n',name,size(prob.A,1),length(prob.sdcone));
    spname = [name 'chk.mat'];
    sdpWrite(prob,spname);
    %Read the problem back from the written files
    probSP = sdpRead(spname);
    %Resolve each
    [~,wSfval(i)] = solve(opti(probSP,opts));
    %Delete intermediate files
    delete(spname);
end

%Display fvals
display([rfval wSfval])
if(norm(rfval-wSfval) > errtol) %LOW tol until Brian corrects repeat problem error
    error('Failed SEDUMI (LP) Sparse Read/Write Check');
end

