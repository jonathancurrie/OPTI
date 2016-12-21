function plotBounds(lb,ub,mplot)
%PLOTBOUNDS Plot bounds on the current figure
%   plotBounds(lb,ub)

%   Copyright (C) 2011 Jonathan Currie (IPL)

xl = xlim; yl = ylim; idx = mplot.idx;
hold on;

%Lower Bounds
if(~isempty(lb))
    if(length(lb) == 1)
        patch([xl(1) max(lb(idx(1)),xl(1)) max(lb(idx(1)),xl(1)) xl(1)],[yl(1) yl(1) yl(2) yl(2)],'y','FaceAlpha',0.3);
    else
        patch([xl(1) xl(2) xl(2) xl(1)],[lb(idx(2)) lb(idx(2)) yl(1) yl(1)],'y','FaceAlpha',0.3)
        patch([lb(idx(1)) lb(idx(1)) xl(1) xl(1)],[yl(1) yl(2) yl(2) yl(1)],'y','FaceAlpha',0.3)
    end
end

%Upper Bounds
if(~isempty(ub))
    if(length(ub) == 1)
        patch([xl(2) min(ub(idx(1)),xl(2)) min(ub(idx(1)),xl(2)) xl(2)],[yl(1) yl(1) yl(2) yl(2)],'y','FaceAlpha',0.3);
    else
        patch([xl(1) xl(2) xl(2) xl(1)],[ub(idx(2)) ub(idx(2)) yl(2) yl(2)],'y','FaceAlpha',0.3)
        patch([ub(idx(1)) ub(idx(1)) xl(2) xl(2)],[yl(1) yl(2) yl(2) yl(1)],'y','FaceAlpha',0.3)
    end
end

hold off;

end

