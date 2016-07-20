function h = plotConfEllipse(thetaOpt,confStats)
%PLOTCONFELLIPSE  Plot Confidence Ellipse/Ellipsoid


scov = confStats.Cov;
nparam = length(confStats.ConfInt);
ndata = size(confStats.ConfBnds.bnds,1);
clim = confStats.Conf;
ci = confStats.ConfInt;

%Find confidence ellipse/ellipsoid
[Ve,D] = eig(scov); e = diag(D); dfe = ndata-nparam; %note update ndata
scale = nparam*(ndata-1)/dfe * rmathlib('qf',clim,nparam,dfe); 

if(length(thetaOpt)==2)
    t = linspace(0,2*pi,5e2)';
    uxc = sqrt(e(1)*scale)*cos(t); uyc = sqrt(e(2)*scale)*sin(t);
    c = ((Ve')\[uxc, uyc]')';
    xc = c(:,1) + thetaOpt(1); yc = c(:,2) + thetaOpt(2);
    modestr = 'Ellipse';
    h(1) = plot(thetaOpt(1),thetaOpt(2),'r*');
    hold on
    h(2) = plot(xc,yc,'b');
    h(3) = plot([-ci(1) -ci(1) ci(1) ci(1) -ci(1)]+thetaOpt(1),[-ci(2) ci(2) ci(2) -ci(2) -ci(2)]+thetaOpt(2),'k');
    hold off
    str = sprintf('%g%%',clim*100);
    legend('Least Squares Optimal Solution',[str ' Confidence ',modestr],['Rectangular ' str ' Confidence Interval'],'location','nw');
else
    optiwarn('opti:plotEllipse','3D Not implemented yet');
%     u = linspace(0,2*pi,100)'; v = linspace(0,pi,length(u))';
%     [U,V] = meshgrid(u,v);    
%     uxc = sqrt(e(1)*scale).*cos(U).*sin(V); uyc = sqrt(e(2)*scale).*sin(U).*cos(V); uzc = sqrt(e(2)*scale)*cos(V);
%     
%     c = ((Ve')\[uxc(:) uyc(:) uzc(:)]')';
%     xc = c(:,1) + thetaOpt(1); yc = c(:,2) + thetaOpt(2); zc = c(:,3) + thetaOpt(3);
% 
%     n = length(u);
%     xc = reshape(xc,n,n);
%     yc = reshape(yc,n,n);
%     zc = reshape(zc,n,n);
% 
%     h(2) = mesh(xc,yc,zc);
    modestr = 'Ellipsoid';
%     hold on
%         h(1) = plot3(thetaOpt(1),thetaOpt(2),thetaOpt(3),'r*');
%     hold off;
end
title(sprintf('Optimal Parameter %g%% Confidence %s',clim*100,modestr));
xlabel('\theta_1'); ylabel('\theta_2'); zlabel('\theta_3');

