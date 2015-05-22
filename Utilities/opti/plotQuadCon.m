function plotQuadCon(Q,l,rl,ru,data)
%PLOTQUADCON Plot Quadratic Constraints on the current figure
%   plotQuadCon(Q,l,r)

%   Copyright (C) 2013 Jonathan Currie (I2C2)

xl = data.xl; yl = data.yl;

hold on;

%Colour
dkg = [0.2 0.2 0.2];

%Determine number of quad constraints
if(iscell(Q))
    no = length(Q);
else
    no = 1;
end
%>= 2014b requires more points...
if(~verLessThan('matlab','8.4'))
    nmul = 2;
else
    nmul = 1;
end

%Generate Constraint Surface Points
[x1,x2] = meshgrid(linspace(xl(1),xl(2),data.npts*nmul),linspace(yl(1),yl(2),data.npts*nmul));

%For each quadratic constraint, plot
for i = 1:no
    %Get Constraint Variables
    if(iscell(Q))
        nQ = Q{i}; nl = l(:,i); nrl = rl(i); nru = ru(i);
    else
        nQ = Q; nl = l; nrl = rl; nru = ru;
    end  
    %Form Constraint Function
    if(data.ndec==1)
        con = @(x) x(1)*nQ*x(1) + nl*x(1);
    else
        con = @(x) x.'*nQ*x + nl.'*x;
    end
    %Plot Each Quad Con as General Nonlinear Constraint
    plotNLCon(con,[],nrl,nru,x1,x2,dkg,data);    
end
hold off;



%OLD CODE
% %Plot Quadratic Inequality Constraints (Inefficient.. need to solve the quadratic)
% [x1,x2] = meshgrid(linspace(xl(1),xl(2),npts),linspace(yl(1),yl(2),npts));
% nox = size(x1);
% noy = size(x2);
% obj = zeros(nox(1),noy(2));
% if(iscell(Q))
%     no = length(Q);
% else
%     no = 1;
% end
% for i = 1:no
%     %get vars
%     if(iscell(Q))
%         nQ = Q{i}; nl = l(:,i); nrl = rl(i); nru = ru(i);
%     else
%         nQ = Q; nl = l; nrl = rl; nru = ru;
%     end    
%     eq = false; sense = 'L';
%     % check for double constraint or equality
%     if(~isinf(nrl) && ~isinf(nru))
%         if(nrl == nru)
%             nr = nru;
%             eq = true;
%         else
%             nr = [nrl nru];
%             sense = 'GL';
%         end
%     elseif(~isinf(nrl))
%         nr = nrl;
%     else
%         nr = nru;
%     end
%     % create surface
%     for n = 1:nox(1)
%         for m = 1:noy(2)
%             x = [x1(n,m) x2(n,m)]';
%             obj(n,m) = x.'*nQ*x + nl.'*x;
%         end
%     end
%     if(eq)
%         contour(x1,x2,obj,'color',[0 0 1],'levellist',nr);
%     else
%         for j = 1:length(nr) %each row constraint (max 2)
%             c = contour(x1,x2,obj,'color',dkg,'levellist',nr(j));
%             %Plot Hatch
%             if(~isempty(c))
%                 %See if we have multiple contours (non-convex or sd)
%                 len = size(c,2)-1;
%                 if(c(2,1) ~= len)
%                     %Build contour array
%                     cstrt = 2; cend = []; n = 2; ind = 1;
%                     while(ind <= len)
%                         ind = ind + c(2,ind) + 1;
%                         cend(n-1) = ind-1; %#ok<AGROW>
%                         cstrt(n) = ind+1;  %#ok<AGROW>
%                         n = n + 1;
%                     end
%                 else
%                     cstrt = 2;
%                     cend = len;
%                 end
%                 %Plot each contour hatch
%                 for n = 1:length(cend)
%                     %Get contour vectors
%                     vecx = diff(c(1,cstrt(n):cend(n)));
%                     vecy = diff(c(2,cstrt(n):cend(n)));
%                     if(isempty(vecx) || isempty(vecy))
%                         continue;
%                     end
%                     %Rotate hatch lines based on infeasible region
%                     xt = [c(1,cstrt(n))+vecy(1) c(2,cstrt(n))-vecx(1)]'; %check rotated -90
%                     fval = xt.'*nQ*xt + nl.'*xt ;      
%                     if(fval <= nr(j) || sense(j) == 'G') %rotate 90
%                         hvecx = -vecy;
%                         hvecy = vecx;
%                     else %rotate -90
%                         hvecx = vecy;
%                         hvecy = -vecx;
%                     end
%                     %Normalize 
%                     av = mean(sqrt(hvecx.^2 + hvecy.^2));
%                     dirs = atan2(hvecy,hvecx);    
%                     hvecx = av*cos(dirs);
%                     hvecy = av*sin(dirs);
%                     %Shift origin
%                     hvecx = c(1,cstrt(n):cend(n)-1) + hvecx;
%                     hvecy = c(2,cstrt(n):cend(n)-1) + hvecy;
%                     %Plot
%                     line([c(1,cstrt(n):cend(n)-1)' hvecx']',[c(2,cstrt(n):cend(n)-1)' hvecy']','Color','k')
%                 end                
%             else
%                 optiwarn('opti:plot','Cannot plot inequality constraint as contour data is empty!');
%             end
%         end
%     end
% end
% 
% hold off;

% CONVEX QC CODE
% %Get contour vectors
% vecx = diff(c(1,2:end));
% vecy = diff(c(2,2:end));
% %Rotate hatch lines based on infeasible region
% xt = [c(1,2)+vecy(1) c(2,2)-vecx(1)]'; %check rotated -90
% fval = xt.'*nQ*xt + nl.'*xt ;   
% if(fval <= nr) %rotate 90
%     hvecx = -vecy;
%     hvecy = vecx;
% else %rotate -90
%     hvecx = vecy;
%     hvecy = -vecx;
% end
% %Normalize
% av = mean(sqrt(hvecx.^2 + hvecy.^2));
% dirs = atan2(hvecy,hvecx);    
% hvecx = av*cos(dirs);
% hvecy = av*sin(dirs);
% %Shift origin
% hvecx = c(1,2:end-1) + hvecx;
% hvecy = c(2,2:end-1) + hvecy;
% %Plot
% line([c(1,2:end-1)' hvecx']',[c(2,2:end-1)' hvecy']','Color','k')