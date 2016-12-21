function plotNonlinCon(prob,data,linear)
%PLOTNONLINCON Plot Nonlinear Constraints on the current figure
%   plotNonlinCon(prob)

%   Copyright (C) 2013 Jonathan Currie (IPL)

if(nargin < 3), linear = false; end
xl = data.xl; yl = data.yl;

hold on;

%Colour
dkg = [0.2 0.2 0.2];
data.linear = linear;

%Ensure we have row format
if(isempty(prob.cl))
    prob = nmix2row(prob);
end

%Determine number of nonlinear constraints
no = length(prob.cl);
%Generate Constraint Surface Points
[x1,x2] = meshgrid(linspace(xl(1),xl(2),data.npts),linspace(yl(1),yl(2),data.npts));
nox = size(x1); noy = size(x2);
%Form Constraint Function
if(prob.sizes.ndec==1)
    con = @(x) prob.nlcon(x(1));
else
    con = prob.nlcon;
end

%Create surface of all constraints
obj = zeros(nox(1),noy(2),no);       
for n = 1:nox(1)
    for m = 1:noy(2)
        x = data.fixval;
        if(data.ndec > 1)
            x(data.idx) = [x1(n,m) x2(n,m)]';
        else
            x(data.idx) = x1(n,m);
        end
        obj(n,m,:) = con(x);
    end
end

%Plot All General Nonlinear Constraints
plotNLCon(con,obj,prob.cl,prob.cu,x1,x2,dkg,data); 

hold off;



% OLD CODE
% if(isempty(prob.nle))
%     prob = nrow2mix(prob,0);
% end
% 
% %Generate Nonlinear Constraint Contours & Plot
% nocon = length(prob.nle);
% if(nocon)
%     [x1,x2] = meshgrid(linspace(xl(1),xl(2),npts),linspace(yl(1),yl(2),npts));
%     nox = size(x1);
%     noy = size(x2);
%     obj = zeros(nox(1),noy(2),nocon);       
%     % create surface of all constraints
%     for n = 1:nox(1)
%         for m = 1:noy(2)
%             x = [x1(n,m) x2(n,m)]';
%             obj(n,m,:) = prob.nlcon(x);
%         end
%     end
%     %Plot them
%     for i = 1:nocon
%         %Inequality
%         if(prob.nle(i)) 
%             c = contour(x1,x2,obj(:,:,i),'color',dkg,'levellist',prob.nlrhs(i));
%             %Plot Hatch
%             if(~isempty(c))
%                 %See if we have multiple contours
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
%                     fval = prob.nlcon(xt);   
%                     if((fval(i) <= prob.nlrhs(i) && prob.nle(i) == -1) || (fval(i) >= prob.nlrhs(i) && prob.nle(i) == 1)) %rotate 90
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
%         else
%             contour(x1,x2,obj(:,:,i),'color',[0 0 1],'levellist',prob.nlrhs(i)); 
%         end         
%     end
% end
% 
% hold off;