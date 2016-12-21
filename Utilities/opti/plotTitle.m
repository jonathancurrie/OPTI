function str = plotTitle(prob,xb,data)
%PLOTTITLE Add Plot Title to Current Axes
    
%   Copyright (C) 2011 Jonathan Currie (IPL)

if(nargin < 2), xb = []; end

%Build Title
str = sprintf('%s Optimization',upper(prob.type));
if(data.ndec==1)
    str = sprintf('%s Function',str); 
else
    str = sprintf('%s Contours',str);
end

if(data.dolog)
    str = sprintf('%s (Log)',str);
end
if(any(strcmpi(prob.type,{'sle','snle','scnle'})) && data.ndec > 1)
    str = sprintf('%s & Equations',str);
end
if(prob.sizes.ncon)
    str = sprintf('%s + Constraints',str);
end
str = sprintf('%s -',str);
if(~isempty(xb) && any(strcmpi(data.mode,{'multi','normal'})))
    if(prob.sense==1)
        str = sprintf('%s Min:',str);
    else
        str = sprintf('%s Max:',str);
    end
    switch(length(xb))
        case 1
            str = sprintf('%s %1.2g',str,xb(1));
        case 2
            str = sprintf('%s [%1.2g; %1.2g]',str,xb(1),xb(2));
        otherwise
            str = [str ' ['];
            for i = 1:length(xb)-1
                str = sprintf('%s %1.2g;',str,xb(i));
            end
            str = sprintf('%s %1.2g]',str,xb(end));
    end                
end
if(~isempty(data.fval) && ~strcmpi(data.mode,'bounded_scale'))
    str = sprintf('%s Fval: %s',str,num2str(data.fval));
else
    str = str(1:end-2);
end

switch(data.mode)
    case 'usex0'
        str = [str ' [At x0]'];
    case {'bounded','bounded_scale'}
        str = [str ' [Within Problem Bounds]'];    
end