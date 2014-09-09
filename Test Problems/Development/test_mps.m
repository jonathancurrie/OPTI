%% TEST MPS READER
clc
tic
prob = coinRead('testLP','mps',1)
toc
% [x,fval,ef,info] = solve(opti(prob,optiset('solver','lpsolve')))

%% TEST MPS WRITER
clc

ftype = 'mps';
opts = optiset('solver','clp');

prob = coinRead('testLP',ftype);
Opt = opti(prob,opts);
[x,fval,ef,stat] = solve(Opt);

tic
coinWrite(prob,'testLPout',ftype)
toc

prob1 = coinRead('testLPout',ftype);
[x1,fval1,ef1,stat1] = solve(opti(prob1,opts));

norm(x-x1)
fval
fval1

delete('testLPout.mps');

%% TEST GMPL READER
% clc
% tic
% prob = coinRead('money','gmpl')
% toc
% % solve(opti(prob))

%% TEST GAMS READER (NOT WORKING - not sure of required GAMS format)
clc
tic
% prob = coinRead('river','gams',1)
toc

%% TEST LP READER
% clc
% tic
% prob = coinRead('prod','lp',1)
% toc

%% P1 - MPS LP
clc
prob = optiReadMPS('testLP');
[x,fval,ef] = solve(opti(prob));
prob = coinRead('testLP','mps');
[x,fval1,ef] = solve(opti(prob));
fval
fval1

%% P2 - MPS LP
clc
prob = optiReadMPS('testLP2')
[x,fval,ef] = solve(opti(prob));

%% P3 - MPS LP
clc
prob = optiReadMPS('testLP3');
[x,fval,ef] = solve(opti(prob));
prob = coinRead('testLP3','mps');
[x,fval1,ef] = solve(opti(prob));
fval
fval1

%% P5 - MPS LP
clc
prob = optiReadMPS('testLP4');
[x,fval,ef] = solve(opti(prob));
prob = coinRead('testLP4','mps');
[x,fval1,ef] = solve(opti(prob));
fval
fval1

%% P6 - MPS MILP
clc
prob = optiReadMPS('testMILP');
[x,fval,ef] = solve(opti(prob));
prob = coinRead('testMILP','mps');
[x,fval1,ef] = solve(opti(prob));
fval
fval1

%% P7 - MPS MILP
clc
prob = optiReadMPS('testMILP2');
[x,fval,ef] = solve(opti(prob));
prob = coinRead('testMILP2','mps');
[x,fval1,ef] = solve(opti(prob));
fval
fval1

%% P8 - QPS QP
clc
prob = optiReadMPS('testQP.qps');
[x,fval,ef] = solve(opti(prob,optiset('solver','ooqp')));
prob = coinRead('testQP.qps');
[x,fval1,ef] = solve(opti(prob,optiset('solver','ooqp')));
fval
fval1

%% P9 - QPS QP
clc
prob = optiReadMPS('testQP2.qps');
[x,fval,ef] = solve(opti(prob,optiset('solver','ooqp')));
prob = coinRead('testQP2','qps');
[x,fval1,ef] = solve(opti(prob,optiset('solver','ooqp')));
fval
fval1

%% Large Read Test
clc
% tic
% prob = optiReadMPS('bienst1')
% toc
tic
prob = coinRead('bienst1.mps')
toc