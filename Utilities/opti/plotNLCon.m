function plotNLCon(fun,con,rl,ru,x1,x2,color,data)
%PLOTNLCON  Plot Hatched Nonlinear Constraint on Current Axes

%   Copyright (C) 2013 Jonathan Currie (IPL)

if(nargin < 6), color = [0.2 0.2 0.2]; end
ss = warning('off','MATLAB:contour:ConstantData');

%Sizes
nox = size(x1); noy = size(x2);
if(isempty(con))
    %Create Surface
    con = zeros(nox(1),noy(2)); x = data.fixval; 
    for n = 1:nox(1)
        for m = 1:noy(2)
            if(data.ndec > 1)
                x(data.idx) = [x1(n,m) x2(n,m)]';
            else
                x(data.idx) = x1(n,m);
            end
            con(n,m) = fun(x);
        end
    end
end

%Determine if we have multiple constraints (user may supply ND con array)
nd = size(con); 
if(length(nd) > 2), ncon = nd(end); else ncon = 1; end

%For each constraint, plot the surface contours
for i = 1:ncon
    %Set Defaults
    eq = false; sense = 'L';
    %Check for double constraint or equality
    if(isempty(rl(i)))
        nr = ru;
    elseif(~isinf(rl(i)) && ~isinf(ru(i)))
        if(rl(i) == ru(i))
            nr = ru;
            eq = true;
        else
            nr = [rl(i) ru(i)];
            sense = 'GL';
        end
    elseif(~isinf(rl(i)))
        nr = rl(i); sense = 'G';
    else
        nr = ru(i);
    end

    %Plot Contours
    %If equality - easy!
    if(eq)
        contour(x1,x2,con(:,:,i),'color',[0 0 1],'levellist',nr(i));
    else %Otherwise inequalities are a bit more work
        for j = 1:length(nr) %each row constraint (max 2)
            if(isfield(data,'linear') && data.linear)
                c = contours(x1,x2,con(:,:,i),[nr(j) nr(j)]); %avoids drawing the contour as well
            else
                c = contour(x1,x2,con(:,:,i),'color',color,'levellist',nr(j));
            end
            %Plot Hatch/Patch
            if(~isempty(c))
                %See if we have multiple contours (non-convex or sd)
                len = size(c,2)-1;
                if(c(2,1) ~= len)
                    %Build contour array
                    cstrt = 2; cend = []; n = 2; ind = 1;
                    while(ind <= len)
                        ind = ind + c(2,ind) + 1;
                        cend(n-1) = ind-1; %#ok<AGROW>
                        cstrt(n) = ind+1;  %#ok<AGROW>
                        n = n + 1;
                    end
                else
                    cstrt = 2;
                    cend = len;
                end
                %Plot each contour hatch/patch
                for n = 1:length(cend)
                    %Get contour vectors
                    xpts = c(1,cstrt(n):cend(n));
                    ypts = c(2,cstrt(n):cend(n));
                    %Ensure no repeated points (otherwise NaNs get through)
                    if(isfield(data,'linear') && data.linear)
                        [~,ix] = unique(xpts); [~,iy] = unique(ypts); ii = intersect(ix,iy);
                        xpts = xpts(ii); ypts = ypts(ii);
                    end
                    vecx = diff(xpts);
                    vecy = diff(ypts);
                    if(isempty(vecx) || isempty(vecy))
                        continue;
                    end
                    %Test point for determining infeasible side
                    xt = data.fixval;
                    if(data.ndec > 1)
                        xt(data.idx) = [c(1,cstrt(n))+vecy(1) c(2,cstrt(n))-vecx(1)]'; %check rotated -90
                    else
                        xt(data.idx) = c(1,cstrt(n))+vecy(1);
                    end
                    fval = fun(xt); fval = fval(i);                                          
                    %Rotate hatch lines based on infeasible region                             
                    if((fval <= nr(j) && sense(j) == 'L') || (fval >= nr(j) && sense(j) == 'G')) %rotate 90
                        hvecx = -vecy;
                        hvecy = vecx;
                    else %rotate -90
                        hvecx = vecy;
                        hvecy = -vecx;
                    end
                    %Normalize 
                    av = mean(sqrt(hvecx.^2 + hvecy.^2));
                    if(all(hvecy==0)), av = 0.4*(av); end
                    dirs = atan2(hvecy,hvecx);                            
                    %If we know is linear, plot as a patch
                    if(isfield(data,'linear') && data.linear)
                        xl = xlim; yl = ylim;
                        %Determine gradient and intercept
                        %vv = vecy./vecx; vv = vv(~isnan(vv) & ~isinf(vv));
                        m = mean(vecy./vecx);
                        ci = mean(ypts-m*xpts);
                        %Find End Points                       
                        xi = (yl - ci) / m;
                        %determine which side to patch
                        md = mean(dirs)*180/pi;
                        if(md >= 90 && md < 180) %top left
                            patch([xi(1) xl(1) xl(1) xi(2)],[yl(1) yl(1) yl(2) yl(2)],'y','FaceAlpha',0.3)                            
                        elseif(md >= 0 && md < 90) %top right
                            patch([xi(1) xl(2) xl(2) xi(2)],[yl(1) yl(1) yl(2) yl(2)],'y','FaceAlpha',0.3)                            
                        elseif(md < 0 && md > -90) % bottom right
                            patch([xi(1) xl(2) xl(2) xi(2)],[yl(1) yl(1) yl(2) yl(2)],'y','FaceAlpha',0.3)
                        else %bottom left
                            patch([xi(1) xl(1) xl(1) xi(2)],[yl(1) yl(1) yl(2) yl(2)],'y','FaceAlpha',0.3)
                        end
                    else %general nonlinear patch
                        hvecx = av*cos(dirs);
                        hvecy = av*sin(dirs);
                        %Shift origin
                        hvecx = c(1,cstrt(n):cend(n)-1) + hvecx;
                        hvecy = c(2,cstrt(n):cend(n)-1) + hvecy;
                        %Plot
                        line([c(1,cstrt(n):cend(n)-1)' hvecx']',[c(2,cstrt(n):cend(n)-1)' hvecy']','Color',color)
                    end
                end                
            else
                optiwarn('OPTI:NonlinearPlot','Cannot plot inequality constraint as contour data is empty!');
            end
        end
    end
end
warning(ss);