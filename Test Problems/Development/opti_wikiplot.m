function opti_wikiplot(mode,data)
% Plots 'Nice' plots for the OPTI Wiki

switch(mode)
    case 'qp'
        qpplot(data.d,data.name);
    case 'riemann'
        riemannPlot;
    case 'wolfram'
        wolframPlot;
    case 'rastrigin'
        rastriginPlot;
        
end


function rastriginPlot
%Objective
fun = @(x) 20 + x(1)^2 + x(2)^2 - 10*(cos(2*pi*x(1)) + cos(2*pi*x(2)));
%Constraints
lb = [-5;-5];
ub = [5;5];

n = 3e2; xmin = -40; xmax = 40;
x = linspace(lb(1),ub(1),n);
y = linspace(lb(2),ub(2),n);
Z = zeros(n,n);
for i = 1:n
    for j = 1:n
        Z(j,i) = fun([x(i),y(j)]);
    end
end
surfc(x,y,Z)
% colormap summer
% shading flat

lightangle(-45,30); shading interp; view(-49,62);
set(gcf,'Renderer','zbuffer');
set(findobj(gca,'type','surface'),...
    'FaceLighting','phong',...
    'AmbientStrength',.3,'DiffuseStrength',.8,...
    'SpecularStrength',.9,'SpecularExponent',25,...
    'BackFaceLighting','unlit')
xlabel('x1'); ylabel('x2'); zlabel('j');

title('Indefinite QP'); 
colormap winter;
xlabel('x1'); ylabel('x2'); title('Rastrigin Function'); 


function wolframPlot
%Objective
obj = @(x) [x(1) - sin(2*x(1) + 3*x(2)) - cos(3*x(1) - 5*x(2));
            x(2) - sin(x(1) - 2*x(2)) + cos(x(1) + 3*x(2))];
fun = @(x) norm(obj(x));
n = 3e2; xmin = -4; xmax = 4;
x = linspace(xmin,xmax,n);
y = linspace(xmin,xmax,n);
Z = zeros(n,n);
for i = 1:n
    for j = 1:n
        Z(j,i) = fun([x(i),y(j)]);
    end
end
surfc(x,y,Z)
lightangle(-45,30); shading interp; view(-49,62);
set(gcf,'Renderer','zbuffer');
set(findobj(gca,'type','surface'),...
    'FaceLighting','phong',...
    'AmbientStrength',.3,'DiffuseStrength',.8,...
    'SpecularStrength',.9,'SpecularExponent',25,...
    'BackFaceLighting','unlit')
xlabel('x1'); ylabel('x2'); zlabel('j');

title('Indefinite QP'); 
colormap winter;
xlabel('x1'); ylabel('x2'); title('Wolfram Global Problem');  



function riemannPlot
n = 1e3;
x = linspace(0.5,2,n); x2 = linspace(0.95,1.05,n);
y = zeros(n,1); y2 = zeros(n,1);
for i = 1:n
    y(i) = RiemannND(x(i));
    y2(i) = RiemannND(x2(i));
end

fig1 = figure(2);
plot(x,y);
title('Riemann Function');
xlabel('x'); ylabel('j');

fig2 = figure(3);
plot(x2,y2);
title('Zoom');

[h1,h2] = inset(fig1,fig2)
set(h2,'xlim',[min(x2) max(x2)])
close(fig1);
close(fig2);




function qpplot(d,name)
n =length(d);
switch(name)
    case 'Positive Semi-Definite QP'
        T = [   -0.9480   -0.5078
   -0.7411   -0.3206];
    case 'Indefinite QP'
        T = [-1.0667    0.3503
    0.9337   -0.0290];
    otherwise
        T=randn(n);
end
[U S V]=svd(T); % no need for symmetrizing
H=U*diag(d)*U';
eig(H);
f = [0;0];
xmin = -1;
xmax = 1;
%Plot      
n = 1e2;
x = linspace(xmin,xmax,n); y = linspace(xmin,xmax,n); Z = zeros(n,n);
for i = 1:n
    for j = 1:n
        xx = [x(i) y(j)]';
        Z(j,i) = 0.5*xx'*H*xx + f'*xx;
    end
end
surfc(x,y,Z);
lightangle(-45,30); shading interp; view(-167,34); 
set(gcf,'Renderer','zbuffer');
set(findobj(gca,'type','surface'),...
    'FaceLighting','phong',...
    'AmbientStrength',.3,'DiffuseStrength',.8,...
    'SpecularStrength',.9,'SpecularExponent',25,...
    'BackFaceLighting','unlit')
xlabel('x1'); ylabel('x2'); zlabel('j');
title(name); 
colormap winter;



function [h_main, h_inset]=inset(main_handle, inset_handle,inset_size)

% The function plotting figure inside figure (main and inset) from 2 existing figures.
% inset_size is the fraction of inset-figure size, default value is 0.35
% The outputs are the axes-handles of both.
% 
% An examle can found in the file: inset_example.m
% 
% Moshe Lindner, August 2010 (C).

if nargin==2
    inset_size=0.35;
end

inset_size=inset_size*.7;
new_fig=figure(1);
main_fig = findobj(main_handle,'Type','axes');
h_main = copyobj(main_fig,new_fig);
set(h_main,'Position',get(main_fig,'Position'))
inset_fig = findobj(inset_handle,'Type','axes');
h_inset = copyobj(inset_fig,new_fig);
ax=get(main_fig,'Position');
set(h_inset,'Position', [.7*ax(1)+ax(3)-inset_size .5*ax(2)+ax(4)-inset_size inset_size inset_size])
