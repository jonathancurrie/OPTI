function plotDataFit(prob,xb,confStats,doLog)
%PLOTDATAFIT Add new Plot with Data + Fit

%   Copyright (C) 2011 Jonathan Currie (I2C2)
    
%Measurement Plot Color & Style
if(~isempty(confStats))
    measC = [0 0 1];
    fitC = [1 0 0];
    measS = 'o';
else
    measC = [0.4 0.4 0.4];
    fitC = [0    0.4470    0.7410];
    measS = 'o';
end

%Use original function
if(isfield(prob.misc,'funorig') && ~isempty(prob.misc.funorig))
    fun = prob.misc.funorig;
else
    fun = prob.fun;
end

%Plot confidence lines if present
if(~isempty(confStats))
    plot(NaN,NaN); %hack to prevent patch removing outer borders
    h = plotConfReg(confStats);
    hold on;
end

if(nargin(fun) == 2)
    if(~isempty(prob.xdata))
        %If xdata is a matrix, then just plot 1:no points
        if(size(prob.xdata,1) > 1 && size(prob.xdata,2) > 1)
            if(size(prob.xdata,1) ~= length(prob.ydata))
                optiwarn('opti:dfit','Cannot plot data fit as size(xdata,1) and length(ydata) are not the same length!');
            end    
            t = 1:length(prob.ydata);
            y = fun(xb,prob.xdata);
            hl(1) = plot(t,prob.ydata,measS,'color',measC); 
            hold on; hl(2) = plot(t,y,'color',fitC); hold off;
        else
            if(length(prob.xdata) ~= length(prob.ydata))
                optiwarn('opti:dfit','Cannot plot data fit as xdata and ydata are not the same length!');
            end    
            %Try plot with Smooth Data
            try
                x = linspace(min(prob.xdata),max(prob.xdata),1e3);
                y = fun(xb,x);
                hl(1) = plot(prob.xdata,prob.ydata,measS,'color',measC);
                hold on; hl(2) = plot(x,y,'color',fitC); hold off;
            catch %just use user supplied data                               
                hl(1) = plot(prob.xdata,prob.ydata,measS,'color',measC);
                hold on; hl(2) = plot(prob.xdata,fun(xb,prob.xdata),'.:'); hold off;
            end                           
        end
        %Add Weights if specified
        if(~isempty(prob.weighting))
            tstr = (['NLS Curve Fit - SSE: ' num2str(sum(((fun(xb,prob.xdata)-prob.ydata).*prob.weighting).^2))]);
        else
            tstr = (['NLS Curve Fit - SSE: ' num2str(sum((fun(xb,prob.xdata)-prob.ydata).^2))]);
        end
        xlabel('x'); ylabel('y');        
    end
else
    t = 1:length(prob.ydata);
    hl(1) = plot(t,prob.ydata,measS,'color',measC); ylabel('y');
    hold on; hl(2) = plot(t,fun(xb),'.:','color',fitC); hold off;
    if(~isempty(prob.weighting))
        tstr = (['NLS Curve Fit - SSE: ' num2str(sum(((fun(xb)-prob.ydata).*prob.weighting).^2))]);
    else
        tstr = (['NLS Curve Fit - SSE: ' num2str(sum((fun(xb)-prob.ydata).^2))]);
    end
    xlabel('x'); 
end

if(~isempty(confStats) && ~isempty(h))
    legend([hl h],'Original Data','NLS Fit',sprintf('%g%% Confidence',confStats.Conf*100));
    title([tstr '  R^2: ' num2str(confStats.Rsquare)]); 
else
    legend('Original Data','NLS Fit');
    title(tstr);
end
