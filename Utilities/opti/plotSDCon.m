function plotSDCon(sdcone,data)
%PLOTSDCON Plot Semidefinite Constraints on the current figure
%   plotSDCon(sdcone)

%   Copyright (C) 2013 Jonathan Currie (IPL)

xl = data.xl; yl = data.yl;
hold on;

%Colour
dkr = [179/255 0.0 70/255];

%Determine number of quad constraints
if(iscell(sdcone))
    no = length(sdcone);
else
    no = 1;
end
%Generate Constraint Surface Points
[x1,x2] = meshgrid(linspace(xl(1),xl(2),data.npts),linspace(yl(1),yl(2),data.npts));

%For each semidefinite constraint, plot
for i = 1:no
    %Get Constraint Variables & Form Constraint Function
    if(iscell(sdcone))
        m = sqrt(size(sdcone{i},1));
        C = reshape(sdcone{i}(:,1),m,m); A0 = reshape(sdcone{i}(:,2),m,m); 
        switch(data.ndec)
            case 1
                con = @(x) min(eig(A0*x(1) - C));
            case 2
                A1 = reshape(sdcone{i}(:,3),m,m);
                con = @(x) min(eig(A0*x(1) + A1*x(2) - C));
            case 3
                A1 = reshape(sdcone{i}(:,3),m,m); A2 = reshape(sdcone{i}(:,4),m,m);
                con = @(x) min(eig(A0*x(1) + A1*x(2) + A2*x(3) - C));
            case 4
                A1 = reshape(sdcone{i}(:,3),m,m); A2 = reshape(sdcone{i}(:,4),m,m); 
                A3 = reshape(sdcone{i}(:,5),m,m);
                con = @(x) min(eig(A0*x(1) + A1*x(2) + A2*x(3) + A3*x(4) - C));
            case 5
                A1 = reshape(sdcone{i}(:,3),m,m); A2 = reshape(sdcone{i}(:,4),m,m); 
                A3 = reshape(sdcone{i}(:,5),m,m); A4 = reshape(sdcone{i}(:,6),m,m);
                con = @(x) min(eig(A0*x(1) + A1*x(2) + A2*x(3) + A3*x(4) + A4*x(5) - C));
            otherwise
                error('Plotting Semidefinite constraints up to 5D is only supported');
        end
    else
        m = sqrt(size(sdcone,1));
        C = reshape(sdcone(:,1),m,m); 
        A0 = reshape(sdcone(:,2),m,m);
        switch(data.ndec)
            case 1
                con = @(x) min(eig(A0*x(1) - C));
            case 2
                A1 = reshape(sdcone(:,3),m,m);
                con = @(x) min(eig(A0*x(1) + A1*x(2) - C));
            case 3
                A1 = reshape(sdcone(:,3),m,m);
                A2 = reshape(sdcone(:,4),m,m);
                con = @(x) min(eig(A0*x(1) + A1*x(2) + A2*x(3) - C));
            case 4
                A1 = reshape(sdcone(:,3),m,m);
                A2 = reshape(sdcone(:,4),m,m);
                A3 = reshape(sdcone(:,5),m,m);
                con = @(x) min(eig(A0*x(1) + A1*x(2) + A2*x(3) + A3*x(4) - C));
            case 5
                A1 = reshape(sdcone(:,3),m,m);
                A2 = reshape(sdcone(:,4),m,m);
                A3 = reshape(sdcone(:,5),m,m);
                A4 = reshape(sdcone(:,6),m,m);
                con = @(x) min(eig(A0*x(1) + A1*x(2) + A2*x(3) + A3*x(4) + A4*x(5) - C));
            otherwise
                error('Plotting Semidefinite constraints up to 5D is only supported');
        end      
    end  
    %Plot Each Quad Con as General Nonlinear Constraint
    plotNLCon(con,[],1e-6,Inf,x1,x2,dkr,data);    
end
hold off;




% OLD CODE 
% %Plot Semidefinite Constraints (Inefficient.. ideas appreciated!)
% [x1,x2] = meshgrid(linspace(xl(1),xl(2),npts),linspace(yl(1),yl(2),npts));
% nox = size(x1);
% noy = size(x2);
% obj = zeros(nox(1),noy(2));
% if(iscell(sdcone))
%     no = length(sdcone);
% else
%     no = 1;
% end
% for i = 1:no
%     %get vars
%     if(iscell(sdcone))
%         C = sdcone{i}(:,1); A0 = sdcone{i}(:,2); A1 = sdcone{i}(:,3);
%     else
%         m = sqrt(size(sdcone,1));
%         C = reshape(sdcone(:,1),m,m); 
%         A0 = reshape(sdcone(:,2),m,m); 
%         A1 = reshape(sdcone(:,3),m,m);
%     end           
%     % create surface
%     for n = 1:nox(1)
%         for m = 1:noy(2)      
%             obj(n,m) = min(eig(A0*x1(n,m) + A1*x2(n,m) - C));
%         end
%     end
%     c = contour(x1,x2,obj,'color',dkr,'levellist',0);
%     %Plot Hatch
%     if(~isempty(c))
%         %See if we have multiple contours (non-convex or sd)
%         len = size(c,2)-1;
%         if(c(2,1) ~= len)
%             %Build contour array
%             cstrt = 2; cend = []; n = 2; ind = 1;
%             while(ind <= len)
%                 ind = ind + c(2,ind) + 1;
%                 cend(n-1) = ind-1; %#ok<AGROW>
%                 cstrt(n) = ind+1;  %#ok<AGROW>
%                 n = n + 1;
%             end
%         else
%             cstrt = 2;
%             cend = len;
%         end
%         %Plot each contour hatch
%         for n = 1:length(cend)
%             %Get contour vectors
%             vecx = diff(c(1,cstrt(n):cend(n)));
%             vecy = diff(c(2,cstrt(n):cend(n)));
%             if(isempty(vecx) || isempty(vecy))
%                 continue;
%             end
%             %Rotate hatch lines based on infeasible region
%             xt = [c(1,cstrt(n))+vecy(1) c(2,cstrt(n))-vecx(1)]'; %check rotated -90
%             fval = all(eig(A0*xt(1) + A1*xt(2) - C) >= 1e-6);      
%             if(fval) %rotate 90
%                 hvecx = -vecy;
%                 hvecy = vecx;
%             else %rotate -90
%                 hvecx = vecy;
%                 hvecy = -vecx;
%             end
%             %Normalize 
%             av = mean(sqrt(hvecx.^2 + hvecy.^2));
%             dirs = atan2(hvecy,hvecx);    
%             hvecx = av*cos(dirs);
%             hvecy = av*sin(dirs);
%             %Shift origin
%             hvecx = c(1,cstrt(n):cend(n)-1) + hvecx;
%             hvecy = c(2,cstrt(n):cend(n)-1) + hvecy;
%             %Plot
%             line([c(1,cstrt(n):cend(n)-1)' hvecx']',[c(2,cstrt(n):cend(n)-1)' hvecy']','Color',dkr)
%         end                
%     else
%         optiwarn('opti:plot','Cannot plot semidefinite constraint as contour data is empty!');
%     end
% end
% 
% hold off;
