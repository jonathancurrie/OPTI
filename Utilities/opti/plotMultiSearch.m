function plotMultiSearch(prob,data)
%PLOTMULTISEARCH  Plot Multi-Start Search Space

lbnd = prob.multi.lbnd;
ubnd = prob.multi.ubnd;
lbnd2 = prob.multi.lbnd2;
ubnd2 = prob.multi.ubnd2;
x_search = prob.multi.x_search;
f_search = prob.multi.f_search;
x_search2 = prob.multi.x_search2;
f_search2 = prob.multi.f_search2;
x_searchx0 = prob.multi.x_searchx0;
f_searchx0 = prob.multi.f_searchx0;

ndec = data.ndec;
idx = data.idx;


%Develop dotsize expression
fs = [f_search;f_search2;f_searchx0]; minfs = abs(min(fs));
% if(minfs < 0), fs = fs + abs(minfs); end
% maxfs = abs(max(fs)); %not sure overall scaling is best?
dotsize = @(fs) 25-fs./(max(fs)./15);

hold on;
%Plot Search 1 Patches
if(ndec == 1)
    plot1DPatch(lbnd,ubnd,0.1);
else
    plotPatch(lbnd(:,idx),ubnd(:,idx),0.1);
end

%Plot Test 1 Points
if(~isempty(x_search) && ~isempty(x_search{1}))
    if(ndec==1)
        plot1DPoint(x_search,prob.objective,dotsize(f_search+minfs),'b.');
    else
        plotPoint(x_search,idx,dotsize(f_search+minfs),'b.');
    end    
end
    
%Plot Search 2 Patches
if(ndec == 1)
    plot1DPatch(lbnd2,ubnd2,0.3);
else
    plotPatch(lbnd2(:,idx),ubnd2(:,idx),0.3);
end

%Plot Test 2 Points
if(~isempty(x_search2) && ~isempty(x_search2{1})) 
    if(ndec==1)
        plot1DPoint(x_search2,prob.objective,dotsize(f_search2+minfs),'m.');
    else
        plotPoint(x_search2,idx,dotsize(f_search2+minfs),'m.');
    end
end

%Plot User Point + Vicinity
if(~isempty(x_searchx0) && ~isempty(x_searchx0{1}))
    if(ndec==1)
        plot1DPoint(x_searchx0,prob.objective,dotsize(f_searchx0+minfs),'g.');
    else
        plotPoint(x_searchx0,idx,dotsize(f_searchx0+minfs),'g.');
    end
end

function plot1DPatch(lbnd,ubnd,alpha)
yl = ylim;
for i = 1:size(lbnd,1)
    patch([lbnd(i) lbnd(i) ubnd(i) ubnd(i)],[yl(1) yl(2) yl(2) yl(1)],zeros(1,4),'facecolor','none','edgealpha',alpha);
end
    
function plotPatch(lbnd,ubnd,alpha)
for i = 1:size(lbnd,1)
    patch([lbnd(i,1) lbnd(i,1) ubnd(i,1) ubnd(i,1)],[lbnd(i,2) ubnd(i,2) ubnd(i,2) lbnd(i,2)],zeros(1,4),'facecolor','none','edgealpha',alpha);
end

function plot1DPoint(xs,obj,ds,marker)
for i = 1:length(xs) 
    plot(xs{i},obj(xs{i}),marker,'markersize',ds(i));
end

function plotPoint(xs,idx,ds,marker)
for i = 1:length(xs) 
    plot(xs{i}(idx(1)),xs{i}(idx(2)),marker,'markersize',ds(i));
end
