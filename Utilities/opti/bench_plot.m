function bench_plot(times,res,solvers,prob)
%BENCH_PLOT  Plot the results of a benchmark test
%
% bench_plot(times,res,solvers,prob) not intended to be called by a user -
% use optiBench.

%   Copyright (C) 2011 Jonathan Currie (Control Engineering)

nS = length(solvers);   % Number of Solvers
nP = length(times{1});  % Number of Test problems

tol = 1e-4;

%Determine maximum time for normalizing
maxT = 0;
for i = 1:nS
    s = sum(times{i});
    if(s > maxT), maxT = s; end
end
% maxT = 1;

[~,srtInd] = sort(times{1});

%Generate Plot Data
X = zeros(2*nP+3,nS); Y = X;
str = cell(nS,1);
for i = 1:nS
    acc = res{i}.acc(srtInd); ind = isinf(acc);
    if(any(ind))
        acc(ind) = 1;
    end
    ind = acc < tol;
    fM = sum(ind)/nP;
    fC = [0;cumsum(ind)/nP;fM]; 
    fT = [0;cumsum(times{i}(srtInd));maxT];
    [X(:,i),Y(:,i)] = stairs(fT,fC);
    str{i} = [upper(solvers{i}) ': ' num2str(fM*100,'%2.2f') '%'];
end

% %Determine maximum time for normalizing
% maxT = 0;
% for i = 1:nS
%     s = sum(times{i});
%     if(s > maxT), maxT = s; end
% end
% 
% %Generate Plot Data
% X = zeros(2*nP+3,nS); Y = X;
% str = cell(nS,1);
% for i = 1:nS
%     acc = res{i}.acc; ind = isinf(acc);
%     if(any(ind))
%         acc(ind) = 1;
%     end
%     ind = acc < tol;
%     fM = sum(ind)/nP;
%     fC = [0;cumsum(ind)/nP;fM]; 
%     fT = [0;cumsum(times{i})/maxT;1];
%     [X(:,i),Y(:,i)] = stairs(fT,fC);
%     str{i} = [upper(solvers{i}) ': ' num2str(fM*100,'%2.2f') '%'];
% end 
    
%Plot
plot(X,Y);
xlabel('Time [s]');
ylabel('Fraction of Problems Solved');
title([upper(prob) ' Benchmark Results for ' num2str(nP) ' Problems']);
axis([0 maxT 0 1.1]);
legend(str,'location','SE');