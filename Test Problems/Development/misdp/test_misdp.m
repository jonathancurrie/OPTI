clc; 
prob = sdpRead('2areas.dat-s');
prob.A = []; prob.b = [];
opts = optiset('display','iter','solver','scip');
Opt = opti(prob,opts)

[~,f,e,i] = solve(Opt)