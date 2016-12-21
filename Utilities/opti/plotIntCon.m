function plotIntCon(prob,data)
%PLOTINTCON Plot Integer Constraints on the current figure
%   plotIntCon(prob)

%   Copyright (C) 2011 Jonathan Currie (IPL)

%Check we have integer variables in this plot
if(all(prob.int.str(data.idx) == 'C'))
    return; %nothing to do!
end

%Colours
intC = [0 114/255 189/255]; %feasible
intIC = [123/255 142/255 147/255]; %infeasible

xl = round(data.xl); yl = round(data.yl);

hold on;

%Get Problem Components
if(~isempty(prob.rl))
    [A,b,Aeq,beq] = row2gen(prob.A,prob.rl,prob.ru);
else
    A = prob.A; b = prob.b;
    Aeq = prob.Aeq; beq = prob.beq;
end
Q = prob.Q; l = prob.l; qrl = prob.qrl; qru = prob.qru;
if(isempty(prob.lb)), lb = -Inf; else lb = prob.lb; end
if(isempty(prob.ub)), ub = Inf; else ub = prob.ub; end

%Check for 1D problem
if(prob.sizes.ndec == 1)
    if(prob.int.ind(1) == -1)
        X = 0:1; Y = [];
    else
        X = (floor(xl(1))+1):(ceil(xl(2))-1); Y = [];
    end
else
    %Generate Grid and check for binary constraints
    if(prob.int.ind(1) == -1)
        x = 0:1;
    else
        try
            x = (floor(xl(1))+1):(ceil(xl(2))-1);
        catch ME
            optiwarn('opti:plotint',['Could not plot integer constraints due to plotting error: ' ME.message]);
            return;
        end
    end
    if(prob.int.ind(2) == -1)
        y = 0:1;
    else
        try
            y = (floor(yl(1))+1):(ceil(yl(2))-1);
        catch ME
            optiwarn('opti:plotint',['Could not plot integer constraints due to plotting error: ' ME.message]);
            return;
        end
    end
    [X,Y] = meshgrid(x,y);
end

%Check for row nonlinear constraints
if(~isempty(prob.cl))
    prob = nrow2mix(prob,0,false);
end

%Work out feasible integer points
idx = zeros(numel(X),1); Xv = X(:); Yv = Y(:);  
for i=1:numel(X)
    if(isempty(Yv))
        xy = Xv(i);
    else
        xy = data.fixval;
        xy(data.idx) = [Xv(i), Yv(i)]';
    end
    if(~isempty(prob.A))
        idx(i) = all(b >= A*xy) && all(xy >= lb) && all(xy <= ub);  
    else
        idx(i) = all(xy >= lb) && all(xy <= ub);  
    end
    if(~isempty(Aeq))
        idx(i) = idx(i) && (abs(Aeq*xy-beq) < 1e-6);
    end
    if(~isempty(Q))
        for n = 1:prob.sizes.nqc
            if(iscell(Q))
                nQ = Q{n}; nl = l(:,n); nrl = qrl(n); nru = qru(n);
            else
                nQ = Q; nl = l; nrl = qrl; nru = qru;
            end
            if(~isinf(nru))
                qu = all(xy'*nQ*xy + nl'*xy <= nru);
            else
                qu = true;
            end
            if(~isinf(nrl))
                ql = all(xy'*-nQ*xy - nl'*xy <= -nrl);
            else
                ql = true;
            end
            idx(i) = idx(i) && ql && qu;
        end
    end
    if(~isempty(prob.nle))
        vals = prob.nlcon(xy);
        leI = find(prob.nle == -1);
        geI = find(prob.nle == 1); 
        eqI = find(prob.nle == 0);
        nlrhs = prob.nlrhs;
        idx(i) = idx(i) && all(vals(leI) <= nlrhs(leI)) && all(vals(geI) >= nlrhs(geI)) && all(abs(vals(eqI)-nlrhs(eqI)) < 1e-6);
    end
end
%Plot points
idx = logical(idx);
if(isempty(Y))
    Y = zeros(size(X));
    for i = 1:length(Y)
        Y(i) = prob.objective(X(i));
    end
end
plot(X(idx),Y(idx),'o','color',intC,'markerface',intC,'markersize',5);
plot(X(~idx),Y(~idx),'o','color',intIC,'markersize',5);

hold off;

end

