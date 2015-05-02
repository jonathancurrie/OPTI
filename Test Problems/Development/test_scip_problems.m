%% SCIP Problem Files
% Worked OK in v3.0.2, no longer in later ones...

%% NLP3
fun = @(x) x(2) + 1e-5*(x(2)-x(1))^2;
lb = [-inf;0];
x0 = [10;1];
fmin = 0;
opts=[];
opts.cipfile = 'nlp3.cip';

opti_scipnl(fun,[],[],[],lb,[],[],[],[],[],x0,opts)

%% NLP26
fun = @(x) (x(1)-x(2))^2 + (x(2)-x(3))^4;
nlcon = @(x) (1+x(2)^2)*x(1) + x(3)^4 - 3;
cl = 0;
cu = 0;
x0 = [-2.6;2;2];
fmin = 0;
opts.cipfile = 'nlp26.cip';

opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);

%% NLP39
fun = @(x) -x(1);
nlcon = @(x) [x(2) - x(1)^3 - x(3)^2;
              x(1)^2 - x(2) - x(4)^2];
cl = [0;0];
cu = [0;0];
x0 = [2;2;2;2];
fmin = -1;
opts = [];
opts.cipfile = 'nlp39.cip';

opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);

%% NLP40
fun = @(x) -x(1)*x(2)*x(3)*x(4);
nlcon = @(x) [x(1)^3 + x(2)^2 - 1;
              x(1)^2*x(4) - x(3);
              x(4)^2 - x(2)];
cl = [0;0;0];
cu = [0;0;0];
x0 = [0.8;0.8;0.8;0.8];
fmin = -0.25;
opts = [];
opts.cipfile = 'nlp40.cip';

opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);

%% NLP47
fun = @(x) (x(1)-x(2))^2 + (x(2)-x(3))^3 + (x(3)-x(4))^4 + (x(4)-x(5))^4;
nlcon = @(x) [x(1) + x(2)^2 + x(3)^3 - 3;
              x(2) - x(3)^2 + x(4) - 1;
              x(1)*x(5) - 1];
cl = [0;0;0];
cu = [0;0;0];
x0 = [2;sqrt(2);-1;2*-sqrt(2);0.5];
fmin = -0.026715;
opts = [];
opts.cipfile = 'nlp40.cip';

opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);