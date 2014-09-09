%% AMPL Testing
% Collection of test problems used for building the AMPL NL interface
clc
if(isempty(which('ampl.exe')))
    haveAmpl = false;
else
    haveAmpl = true;
end

%% LP
clc
prob = amplRead('diet.nl');

Opt = opti(prob,optiset('solver','auto'))
[x,fval] = solve(Opt)


%% LP2
clc
prob = amplRead('prod.nl');

Opt = opti(prob,optiset('solver','clp','display','iter'))
[x,fval,ef,info] = solve(Opt)

%% LP3
clc
prob = amplRead('oil.nl');

Opt = opti(prob)
[x,fval,ef,info] = solve(Opt)

%% LP4
clc
prob = amplRead('mod8.nl');

Opt = opti(prob)
[x,fval,ef,info] = solve(Opt)

%% LP4
clc
prob = amplRead('train1.nl');

Opt = opti(prob)
[x,fval,ef,info] = solve(Opt)

%% MILP
clc
prob = amplRead('multi.nl');

Opt = opti(prob)
[x,fval,ef,info] = solve(Opt)

%% MILP1
clc
prob = amplRead('multmip1.nl');

Opt = opti(prob)
[x,fval,ef,info] = solve(Opt)

%% MILP2
clc
prob = amplRead('multmip2.nl');

Opt = opti(prob)
[x,fval,ef,info] = solve(Opt)

%% MILP3
clc
prob = amplRead('multmip3.nl');

Opt = opti(prob)
[x,fval,ef,info] = solve(Opt)

%% QP
clc
if(haveAmpl)
    prob = amplRead('testQP.mod');
end

Opt = opti(prob)
[x,fval,ef,info] = solve(Opt)

%% QCQP
clc
if(haveAmpl)
    prob = amplRead('testQCQP.mod',[],[],1);
    Opt = opti(prob,optiset('solver','ipopt'))
    [x,fval,ef,info] = solve(Opt)
end

%% QCQP2
clc
if(haveAmpl)
    prob = amplRead('testQCQP2.mod');
    Opt = opti(prob,optiset('solver','ipopt'))
    [x,fval,ef,info] = solve(Opt)
end

%% QCQP5
clc
% clear all
prob = amplRead('hs116.nl',[],[],1); %non-convex QC so solve as NLP

Opt = opti(prob)
[x,fval,ef,stat] = solve(Opt)

%% MIQP
clc
prob = amplRead('testMIQP.nl');

Opt = opti(prob)
[x,fval,ef,info] = solve(Opt)

%% MIQP2
clc
prob = amplRead('testMIQP2.nl');

Opt = opti(prob)
[x,fval,ef,info] = solve(Opt)

%% MIQP3
clc
prob = amplRead('testMIQP3.nl');

Opt = opti(prob)
[x,fval,ef,info] = solve(Opt)
plot(Opt,4)

h = gca;
xl = get(h,'XLabel'); yl = get(h,'YLabel'); tl = get(h,'Title');
set(xl,'fontsize',12,'fontweight','bold')
set(yl,'fontsize',12,'fontweight','bold')
set(tl,'fontsize',12,'fontweight','bold')


%% NLP
clc
prob = amplRead('hs100.nl');

Opt = opti(prob)
[x,fval] = solve(Opt)

asl('close')

%% NLP2
clc
prob = amplRead('ch3.nl');

Opt = opti(prob)
[x,fval] = solve(Opt)

asl('close')

%% NLP3
clc
prob = amplRead('p2gon.nl');

Opt = opti(prob)
[x,fval,ef,stat] = solve(Opt)

asl('close')

%% NLP4 (bounded)
clc
prob = amplRead('camel1.nl');

Opt = opti(prob)
[x,fval,ef,stat] = solve(Opt)

asl('close')

%% NLP5 (bounded)
clc
prob = amplRead('griewank.nl');

Opt = opti(prob)
[x,fval,ef,stat] = solve(Opt)

asl('close')

%% NLP6 
clc
if(haveAmpl)
    prob = amplRead('hs71.mod');    
    Opt = opti(prob)
    [x,fval,ef,stat] = solve(Opt)
    asl('close')
end

%% NLP6 (lin con)
clc
if(haveAmpl)
    prob = amplRead('hs71lin.mod'); 
    Opt = opti(prob)
    [x,fval,ef,stat] = solve(Opt)
    asl('close')
end

%% NLP7 (uncon)
clc
prob = amplRead('camel2.nl');

Opt = opti(prob)
[x,fval,ef,stat] = solve(Opt)

asl('close')

%% MINLP
clc
prob = amplRead('weapon.nl');
opts = optiset('display','iter','solverOpts',bonminset('algorithm','B-BB'));

Opt = opti(prob,opts)
[x,fval,ef,stat] = solve(Opt)

asl('close')

%% QP constant
if(haveAmpl)
    prob = amplRead('testQP.mod')
    prob2 = amplRead('testQPcon.mod')
    opts = optiset('display','iter');
    Opt = opti(prob,opts)
    [x,fval,ef,stat] = solve(Opt)

    Opt2 = opti(prob2,opts)
    [x,fval,ef,stat] = solve(Opt2)
end

