classdef optifit
%OPTIFIT  OPTI Object for Data Fitting / Model Regression / Parameter Optimization.
%
%   ofit = optifit(xdata,ydata,model,x0,lb,ub,wts,opts)
%
%   ofit = optifit(xdata,ydata,model) accepts curve fitting data xdata and
%   ydata, and attempts to fit the model parameters (theta) to the data. 
%   Valid model structures are shown below.
%
%   ofit = optifit(xdata,{ydata zdata},model) accepts surface fitting data
%   xdata, ydata and zdata and attempts to fit the model parameters to the
%   data. Valid model structures are shown below.
%
%   ofit = optifit(xdata,ydata,model,x0) allows the user to specify a
%   general nonlinear function as the model, together with an initial guess 
%   x0. This calling form works for both curve and surface fitting.
%
%   ofit = optifit(xdata,...,x0,lb,ub) allows the user to specify bounds on
%   the parameters. All bounds should be finite, if specified.
%
%   ofit = optifit(xdata,...,ub,wts) allows the user to specify a vector of
%   weights for each data point you are fitting to. The default weight is
%   1.
%
%   ofit = optifit(xdata,...,wts,opts) allows the user to specify optiset
%   options to pass to the optimizer.
%
%
%   Valid Model Structures:
%   - Curve Fitting
%       - poly1:  a*x + b
%       - poly2:  a*x^2 + b*x + c
%       - poly3:  a*x^3 + b*x^2 + c*x + d
%       - poly4:  a*x^4 + b*x^3 + c*x^2 + d*x + e
%       - poly5:  a*x^5 + b*x^4 + c*x^3 + d*x^2 + e*x + f
%       - power1: a*x^b
%       - power2: a*x^b + c
%       - exp1:   a*exp(b*x)
%       - exp2:   a*exp(b*x) + c
%       - exp3:   a*exp(b*x) + c*exp(d*x)
%       - rat01:  a / (x + b)
%       - rat02:  a / (x^2 + b*x + c)
%       - rat03:  a / (x^3 + b*x^2 + c*x + d)
%       - rat11:  (a*x + b) / (x + c)
%       - rat12:  (a*x + b) / (x^2 + c*x + d)
%       - rat13:  (a*x + b) / (x^3 + c*x^2 + d*x + e)
%       - rat21:  (a*x^2 + b*x + c) / (x + d)
%       - rat22:  (a*x^2 + b*x + c) / (x^2 + d*x + e)
%       - rat23:  (a*x^2 + b*x + c) / (x^3 + d*x^2 + e*x + f)
%       - sin1:   a*sin(b*x + c) + d
%       - sin2:   a*sin(b*x + c) + d*sin(e*x + f) + g
%       - sin3:   a*sin(b*x + c) + d*sin(e*x + f) + g*sin(h*x + i) + j
%       - auto:   Will attempt to find the best structure from the above (up to 2nd order).
%
%   - Surface Fitting
%       - poly11: a*x + b*y + c
%       - poly12: a*y^2 + b*x*y + c*x + d*y + e
%       - poly13: a*y^3 + b*y^2 + c*x*y^2 + d*x*y + e*x + f*y + g
%       - poly21: a*x^2 + b*x*y + c*x + d*y + e
%       - poly22: a*x^2 + b*y^2 + c*x*y + d*x + e*y + f
%       - poly23: a*y^3 + b*x^2 + c*y^2 + d*x^2*y + e*x*y^2 + f*x*y + g*x + h*y + i
%       - poly31: a*x^3 + b*x^2 + c*x^2*y + d*x*y + e*x + f*y + g
%       - poly32: a*x^3 + b*x^2 + c*y^2 + d*x^2*y + e*x*y^2 + f*x*y + g*x + h*y + i
%       - poly33: a*x^3 + b*y^3 + c*x^2 + d*y^2 + e*x^2*y + f*x*y^2 + g*x*y + h*x + i*y + j
%
%
%   See also opti opti.solve
%
%   Copyright (C) 2011-2014 Jonathan Currie (www.inverseproblem.co.nz)
    
    
    properties (SetAccess = private)
        theta
        stats
        data
        fun
        grad
        optiObj
        is3D        
    end
    
    methods
        
        %-- Constructor & Solver --%
        function ofit = optifit(xdata,ydata,model,x0,lb,ub,wts,opts)
            %OPTIFIT Constructor
            if(nargin < 8), opts = optiset; end
            if(nargin < 7), wts = []; end
            if(nargin < 6), ub = []; end
            if(nargin < 5), lb = []; end
            if(nargin < 4), x0 = []; end
            if(nargin < 3 || isempty(model)), model = 'auto'; end
            %Pre-process data
            [xd,yd,zd,wts,ofit.data] = optifit.ppData(xdata,ydata,wts);
            %Curve Fit
            if(isempty(zd))
                ofit.is3D = false;
                %If auto model, attempt to find best structure
                if(ischar(model) && strcmpi(model,'auto'))
                    [ofit.theta,ofit.optiObj,ofit.stats,ofit.fun,ofit.grad] = optifit.guessModelandFit(xd,yd,x0,lb,ub,wts,opts);
                else
                    %Fit Model
                    [ofit.theta,ofit.optiObj,ofit.stats,ofit.fun,ofit.grad] = optifit.modelFit(xd,yd,model,x0,lb,ub,wts,opts);
                end
            %Surface Fit
            else
                ofit.is3D = true;
                %If auto model, attempt to find best structure
                if(ischar(model) && strcmpi(model,'auto'))
                    error('auto structure ID not implemented for 3D');
                else
                    %Fit Model
                    [ofit.theta,ofit.optiObj,ofit.stats,ofit.fun,ofit.grad] = optifit.model3DFit(xd,yd,zd,model,x0,lb,ub,wts,opts);
                end              
            end            
        end
        
        %-- Display --%
        function display(ofit)
            optiPrintFitStats(ofit.stats,ofit.theta);
        end
        
        %-- Plot --%
        function plot(ofit)
            %PLOT  Plots the fitted curve/surface 
            if(ofit.is3D)
                plot3(ofit);
            else
                plot(ofit.optiObj);
                title(['Curve Fit - SSE: ' num2str(ofit.stats.SSE) ' R^2: ' num2str(ofit.stats.Rsquare)]); 
            end
            axis('tight');
        end
        
        %-- 3D Plot --%
        function plot3(ofit)
            if(~ofit.is3D)
                error('This function is only for plotting 3D surface fits');
            end     
            %Restore original data shape, including NaNs but without Infs
            xd = NaN(size(ofit.data.idx)); yd = xd; zd = yd; zm = yd;
            xd(ofit.data.idx) = ofit.data.xd;
            yd(ofit.data.idx) = ofit.data.yd;
            zd(ofit.data.idx) = ofit.data.zd;            
            zm(ofit.data.idx) = ofit.optiObj.prob.fun(ofit.theta);
            if(size(zm,1) > 1 && size(zm,2) > 1)
                surf(xd,yd,zm);
                colormap('pink');
            else
                plot3(xd,yd,zm,'*');
            end
            hold on;
            plot3(xd,yd,zd,'k.','markersize',15)
            hold off;
            xlabel('xdata'); ylabel('ydata'); zlabel('zdata');            
            title(['Surface Fit - SSE: ' num2str(ofit.stats.SSE) ' R^2: ' num2str(ofit.stats.Rsquare)]); 
        end
        
        %-- Function Evaluation --%
        function f = feval(ofit,theta,xdata,ydata)
            %FEVAL  Evaluate Model Function with Specified Parameters and
            %optionally xdata/ydata
            %
            %   f = feval(ofit,theta,xdata,ydata)
            
            switch(nargin(ofit.fun))
                case 1
                    f = ofit.fun(theta);
                case 2
                    if(nargin < 3)
                        f = ofit.fun(theta,ofit.data.xd);
                    else
                        f = ofit.fun(theta,xdata);
                    end
                case 3
                    if(nargin < 3)
                        f = ofit.fun(theta,ofit.data.xd,ofit.data.yd);
                    elseif(nargin < 4)
                        f = ofit.fun(theta,xdata,ofit.data.yd);
                    else
                        f = ofit.fun(theta,xdata,ydata);
                    end
            end            
        end
            
        %-- Gradient Evaluation --%
        function g = geval(ofit,theta,xdata,ydata)
            %GEVAL  Evaluate Model Gradient with Specified Parameters and
            %optionally xdata/ydata
            %
            %   g = geval(ofit,theta,xdata,ydata)
            
            switch(nargin(ofit.grad))
                case 1
                    g = ofit.grad(theta);
                case 2
                    if(nargin < 3)
                        g = ofit.grad(theta,ofit.data.xd);
                    else
                        g = ofit.grad(theta,xdata);
                    end
                case 3
                    if(nargin < 3)
                        g = ofit.grad(theta,ofit.data.xd,ofit.data.yd);
                    elseif(nargin < 4)
                        g = ofit.grad(theta,xdata,ofit.data.yd);
                    else
                        g = ofit.grad(theta,xdata,ydata);
                    end
            end
        end
        
        %-- Print Model to Command Line --%
        function str = print(ofit)
            %PRINT Print Model to Command Line [Unfinished]
            %
            %   print(ofit)
            
            str = ofit.stats.ModelStructure;
            for i = length(ofit.theta):-1:1
                t = ['t' int2str(i)];
                str = strrep(str,t,num2str(ofit.theta(i),'%1.12g'));
            end
            idx = strfind(str,'[');
            str = str(1:idx(1)-1);
            if(~nargout)
                disp(str);
            end
        end
            
        
    end
    
    methods (Static, Access=private)
        
        %-- Pre-process Data --%
        function [xd,yd,zd,wts,data] = ppData(xdata,ydata,wts)
           
            %See if we have 3D data
            if(iscell(ydata))
                if(length(ydata) ~= 2)
                    error('For 3D fitting, it is expected that ydata contains {ydata,zdata}');
                end
                xd = xdata; yd = ydata{1}; zd = ydata{2};
                xvec=0;yvec=0;zvec=0;
                %Check for vectors
                if(size(xd,1)==1 || size(xd,2)==1)
                    if(size(xd,2) > 1), xd = xd'; end
                    xvec=1;
                end
                if(size(yd,1)==1 || size(yd,2)==1)
                    if(size(yd,2) > 1), yd = yd'; end
                    yvec=1;
                end
                if(size(zd,1)==1 || size(zd,2)==1)
                    if(size(zd,2) > 1), zd = zd'; end
                    zvec=1;
                end
                if(xor(xvec,yvec)), error('Both xdata and ydata must be vectors, if one is'); end
                if(zvec && (~xvec || ~yvec)), error('zdata must be a matrix if xdata or ydata is a matrix'); end
                if(xvec && yvec && ~zvec)
                    %if not the same number of elements, try meshgrid
                    if(numel(xd) ~= numel(zd))
                        [X,Y] = meshgrid(xd,yd);
                        if(numel(X) ~= numel(zd))
                            error('Dimensions of inputs don''t match! Tried [X,Y] = meshgrid(xdata,ydata), but numel(X) ~= numel(zdata)!');
                        end
                        xd = X; yd = Y;
                    end
                end                                    
                
                if(numel(xd) ~= numel(yd)), error('The number of elements in xdata and ydata don''t match'); end
                if(numel(xd) ~= numel(zd)), error('The number of elements in xdata and zdata don''t match'); end
                if(issparse(xd) || issparse(yd) || ~isreal(xd) || ~isreal(yd) || ~isreal(zd) || ~isreal(zd))
                    error('Fitting data must be real, dense vectors/matrices');
                end
                %Remove NaNs, Infs, returns as columns
                idx = isnan(xd) | isnan(yd) | isnan(zd) | isinf(xd) | isnan(yd) | isinf(zd);
                xd = xd(~idx);
                yd = yd(~idx);
                zd = zd(~idx);
                data = struct('xd',xd,'yd',yd,'zd',zd,'idx',~idx);

            else %2D data            
                %Check lengths, data type
                if(numel(xdata) ~= numel(ydata)), error('The length of xdata must equal ydata'); end                
                xd = xdata(:); yd = ydata(:);
                if(issparse(xd) || issparse(yd) || ~isreal(xd) || ~isreal(yd))
                    error('Fitting data must be real, dense vectors');
                end
                %Remove NaNs, Infs
                idx = isnan(xd) | isnan(yd) | isinf(xd) | isnan(yd);
                xd = xd(~idx);
                yd = yd(~idx);                
                zd = []; 
                data = struct('xd',xd,'yd',yd,'zd',[],'idx',~idx);
            end
            
            %Process weights
            if(~isempty(wts))                
                if(any(isnan(wts)) || any(isinf(wts)) || issparse(wts) || ~isreal(wts))
                    error('Weights must be finite, not NaN, dense, real values');
                end
                data.wts = wts;
                wts = wts(~idx);
                if(numel(wts) ~= numel(yd))
                    error('You must supply the same number of weights as fitting values');
                end                
            else
                data.wts = [];
            end
        end
        
        %-- Fit 2D Model --%
        function [theta,Opt,stats,fun,grad] = modelFit(xdata,ydata,model,x0,lb,ub,wts,opts)
           
            stats = [];
            %Standard Models - Use function knowledge to generate reasonable initial guess.
            if(ischar(model))
                s = warning('off','MATLAB:rankDeficientMatrix');
                s1 = warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
                switch(lower(model))                    
                    case 'poly1'
                        theta0 = polyfit(xdata,ydata,1); 
                        fun = @(theta,xdata) theta(1)*xdata + theta(2);
                        grad = @(theta,xdata) [xdata ones(size(xdata))];
                        str = 't1*xd + t2';
                        
                    case 'poly2'
                        theta0 = polyfit(xdata,ydata,2);
                        fun = @(theta,xdata) theta(1)*xdata.^2 + theta(2)*xdata + theta(3);
                        grad = @(theta,xdata) [xdata.^2 xdata ones(size(xdata))];
                        str = 't1*xd^2 + t2*xd + t3';
                        
                    case 'poly3'
                        theta0 = polyfit(xdata,ydata,3);  
                        fun = @(theta,xdata) theta(1)*xdata.^3 + theta(2)*xdata.^2 + theta(3)*xdata + theta(4);
                        grad = @(theta,xdata) [xdata.^3 xdata.^2 xdata ones(size(xdata))];
                        str = 't1*xd^3 + t2*xd^2 + t3*xd + t4';
                        
                    case 'poly4'
                        theta0 = polyfit(xdata,ydata,4);                        
                        fun = @(theta,xdata) theta(1)*xdata.^4 + theta(2)*xdata.^3 + theta(3)*xdata.^2 + theta(4)*xdata + theta(5);
                        grad = @(theta,xdata) [xdata.^4 xdata.^3 xdata.^2 xdata ones(size(xdata))];
                        str = 't1*xd^4 + t2*xd^3 + t3*xd^2 + t4*xd + t5';
                        
                    case 'poly5'
                        theta0 = polyfit(xdata,ydata,5);                        
                        fun = @(theta,xdata) theta(1)*xdata.^5 + theta(2)*xdata.^4 + theta(3)*xdata.^3 + theta(4)*xdata.^2 + theta(5)*xdata + theta(6);
                        grad = @(theta,xdata) [xdata.^5 xdata.^4 xdata.^3 xdata.^2 xdata ones(size(xdata))];
                        str = 't1*xd^5 + t2*xd^4 + t3*xd^3 + t4*xd^2 + t5*xd + t6';
                        
                    case 'power1' %a*x^b
                        if(any(xdata <= 0)), error('A power model cannot be fitted to data <= 0'); end
                        if(isempty(x0)), theta0 = optifit.genTheta0(model,xdata,ydata,wts); else theta0 = x0; end
                        fun = @(theta,xdata) theta(1)*xdata.^theta(2);
                        grad = @(theta,xdata) [xdata.^theta(2) theta(1)*xdata.^theta(2) .* log(xdata)];
                        str = 't1*xd^t2';
                        
                    case 'power2' %a*x^b + c
                        if(any(xdata <= 0)), error('A power model cannot be fitted to data <= 0'); end
                        if(isempty(x0)), theta0 = optifit.genTheta0(model,xdata,ydata,wts); else theta0 = x0; end
                        fun = @(theta,xdata) theta(1)*xdata.^theta(2) + theta(3);
                        grad = @(theta,xdata) [xdata.^theta(2) theta(1)*xdata.^theta(2) .* log(xdata) ones(size(xdata))];
                        str = 't1*xd^t2 + t3';
                        
                    case 'exp1' %a*exp(b*x)
                        if(isempty(x0)), theta0 = optifit.genTheta0(model,xdata,ydata,wts); else theta0 = x0; end
                        fun = @(theta,xdata) theta(1)*exp(theta(2)*xdata);
                        grad = @(theta,xdata) [exp(theta(2)*xdata) theta(1)*xdata.*exp(theta(2)*xdata)];
                        str = 't1*exp(xd*t2)';
                        
                    case 'exp2' %a*exp(b*x) + c
                        if(isempty(x0)), theta0 = optifit.genTheta0(model,xdata,ydata,wts); else theta0 = x0; end
                        fun = @(theta,xdata) theta(1)*exp(theta(2)*xdata) + theta(3);
                        grad = @(theta,xdata) [exp(theta(2)*xdata) theta(1)*xdata.*exp(theta(2)*xdata) ones(size(xdata))];                            
                        str = 't1*exp(xd*t2) + t3';
                        
                    case 'exp3' %a*exp(b*x) + c*exp(d*x)
                        if(isempty(x0)), theta0 = optifit.genTheta0(model,xdata,ydata,wts); else theta0 = x0; end
                        fun = @(theta,xdata) theta(1)*exp(theta(2)*xdata) + theta(3)*exp(theta(4)*xdata);
                        grad = @(theta,xdata) [exp(theta(2)*xdata) theta(1)*xdata.*exp(theta(2)*xdata) exp(theta(4)*xdata) theta(3)*xdata.*exp(theta(4)*xdata)];
                        str = 't1*exp(xd*t2) + t3*exp(xd*t4)';  
                        
                    case 'rat01' %a / (x + b)
                        theta0 = [ones(size(ydata)) -ydata]\(ydata.*xdata);
                        fun = @(theta,xdata) theta(1) ./ (xdata + theta(2));
                        grad = @(theta,xdata) [1./(theta(1) + xdata) -theta(1)./(theta(2) + xdata).^2];
                        str = 't1 / (x + t2)';
                        
                    case 'rat02' %a / (x^2 + b*x + c)
                        theta0 = [ones(size(ydata)) -ydata.*xdata -ydata]\(ydata.*xdata.^2);
                        fun = @(theta,xdata) theta(1) ./ (xdata.^2 + theta(2)*xdata + theta(3));
                        grad = @(theta,xdata) [1./(xdata.^2 + theta(2)*xdata + theta(3)), -(theta(1)*xdata)./(xdata.^2 + theta(2)*xdata + theta(3)).^2, -theta(1)./(xdata.^2 + theta(2)*xdata + theta(3)).^2];
                        str = 't1 / (x^2 + t2*x + t3)';
                        
                    case 'rat03' %a / (x^3 + b*x^2 + c*x + d)
                        theta0 = [ones(size(ydata)) -ydata.*xdata.^2 -ydata.*xdata -ydata]\(ydata.*xdata.^3);
                        fun = @(theta,xdata) theta(1) ./ (xdata.^3 + theta(2)*xdata.^2 + theta(3)*xdata + theta(4));
                        grad = @(theta,xdata) [1./(xdata.^3 + theta(2)*xdata.^2 + theta(3)*xdata + theta(4)), -(theta(1)*xdata.^2)./(xdata.^3 + theta(2)*xdata.^2 + theta(3)*xdata + theta(4)).^2,...
                                               (-theta(1)*xdata)./(xdata.^3 + theta(2)*xdata.^2 + theta(3)*xdata + theta(4)).^2,-theta(1)./(xdata.^3 + theta(2)*xdata.^2 + theta(3)*xdata + theta(4)).^2];
                        str = 't1 / (x^3 + t2*x^2 + t3*x + t4)';
                        
                    case 'rat11' %(a*x + b) / (x + c)
                        theta0 = [xdata ones(size(xdata)) -ydata]\(ydata.*xdata);
                        fun = @(theta,xdata) (theta(1)*xdata + theta(2)) ./ (xdata + theta(3));
                        grad = @(theta,xdata) [xdata./(theta(3) + xdata), 1./(theta(3) + xdata), -(theta(2) + theta(1)*xdata)./(theta(3) + xdata).^2];
                        str = '(t1*x + t2) / (x + t3)';
                        
                    case 'rat12' %(a*x + b) / (x^2 + c*x + d)
                        theta0 = [xdata ones(size(xdata)) -ydata.*xdata -ydata]\(ydata.*xdata.^2);
                        fun = @(theta,xdata) (theta(1)*xdata + theta(2)) ./ (xdata.^2 + theta(3)*xdata + theta(4));
                        grad = @(theta,xdata) [xdata./(xdata.^2 + theta(3)*xdata + theta(4)), 1./(xdata.^2 + theta(3)*xdata + theta(4)),...
                                               -(xdata.*(theta(2) + theta(1)*xdata))./(xdata.^2 + theta(3)*xdata + theta(4)).^2,...
                                               -(theta(2) + theta(1)*xdata)./(xdata.^2 + theta(3)*xdata + theta(4)).^2];
                        str = '(t1*x + t2) / (x^2 + t3*x + t4)';
                        
                    case 'rat13' %(a*x + b) / (x^3 + c*x^2 + d*x + e)
                        theta0 = [xdata ones(size(xdata)) -ydata.*xdata.^2 -ydata.*xdata -ydata]\(ydata.*xdata.^3);
                        fun = @(theta,xdata) (theta(1)*xdata + theta(2)) ./ (xdata.^3 + theta(3)*xdata.^2 + theta(4)*xdata + theta(5));
                        grad = @(theta,xdata) [xdata./(xdata.^3 + theta(3)*xdata.^2 + theta(4)*xdata + theta(5)), 1./(xdata.^3 + theta(3)*xdata.^2 + theta(4)*xdata + theta(5)),...
                                               -(xdata.^2.*(theta(2) + theta(1)*xdata))./(xdata.^3 + theta(3)*xdata.^2 + theta(4)*xdata + theta(5)).^2,...
                                               -(xdata.*(theta(2) + theta(1)*xdata))./(xdata.^3 + theta(3)*xdata.^2 + theta(4)*xdata + theta(5)).^2,...
                                               -(theta(2) + theta(1)*xdata)./(xdata.^3 + theta(3)*xdata.^2 + theta(4)*xdata + theta(5)).^2];
                        str = '(t1*x + t2) / (x^3 + t3*x^2 + t4*x + t5)';
                        
                    case 'rat21' %(a*x^2 + b*x + c) / (x + d)
                        theta0 = [xdata.^2 xdata ones(size(xdata)) -ydata]\(ydata.*xdata);
                        fun = @(theta,xdata) (theta(1)*xdata.^2 + theta(2)*xdata + theta(3)) ./ (xdata + theta(4));
                        grad = @(theta,xdata) [xdata.^2./(theta(4) + xdata), xdata./(theta(4) + xdata), 1./(theta(4) + xdata), -(theta(1)*xdata.^2 + theta(2)*xdata + theta(3))./(theta(4) + xdata).^2];
                        str = '(t1*x^2 + t2*x + t3) / (x + t4)';
                        
                    case 'rat22' %(a*x^2 + b*x + c) / (x^2 + d*x + e)
                        theta0 = [xdata.^2 xdata ones(size(xdata)) -ydata.*xdata -ydata]\(ydata.*xdata.^2);
                        fun = @(theta,xdata) (theta(1)*xdata.^2 + theta(2)*xdata + theta(3)) ./ (xdata.^2 + theta(4)*xdata + theta(5));
                        grad = @(theta,xdata) [xdata.^2./(xdata.^2 + theta(4)*xdata + theta(5)), xdata./(xdata.^2 + theta(4)*xdata + theta(5)), 1./(xdata.^2 + theta(4)*xdata + theta(5)),...
                                               -(xdata.*(theta(1)*xdata.^2 + theta(2)*xdata + theta(3)))./(xdata.^2 + theta(4)*xdata + theta(5)).^2,...
                                               -(theta(1)*xdata.^2 + theta(2)*xdata + theta(3))./(xdata.^2 + theta(4)*xdata + theta(5)).^2];
                        str = '(t1*x^2 + t2*x + t3) / (x^2 + t4*x + t5)';
                        
                    case 'rat23' %(a*x^2 + b*x + c) / (x^3 + d*x^2 + e*x + f)
                        theta0 = [xdata.^2 xdata ones(size(xdata)) -ydata.*xdata.^2 -ydata.*xdata -ydata]\(ydata.*xdata.^3);
                        fun = @(theta,xdata) (theta(1)*xdata.^2 + theta(2)*xdata + theta(3)) ./ (xdata.^3 + theta(4)*xdata.^2 + theta(5)*xdata + theta(6));
                        grad = @(theta,xdata) [xdata.^2./(xdata.^3 + theta(4)*xdata.^2 + theta(5)*xdata + theta(6)), xdata./(xdata.^3 + theta(4)*xdata.^2 + theta(5)*xdata + theta(6)), ...
                                               1./(xdata.^3 + theta(4)*xdata.^2 + theta(5)*xdata + theta(6)),...
                                               -(xdata.^2.*(theta(1)*xdata.^2 + theta(2)*xdata + theta(3)))./(xdata.^3 + theta(4)*xdata.^2 + theta(5)*xdata + theta(6)).^2,...
                                               -(xdata.*(theta(1)*xdata.^2 + theta(2)*xdata + theta(3)))./(xdata.^3 + theta(4)*xdata.^2 + theta(5)*xdata + theta(6)).^2,...
                                               -(theta(1)*xdata.^2 + theta(2)*xdata + theta(3))./(xdata.^3 + theta(4)*xdata.^2 + theta(5)*xdata + theta(6)).^2];
                        str = '(t1*x^2 + t2*x + t3) / (x^3 + t4*x^2 + t5*x + t6)';
                        
                    case 'sin1' %a*sin(b*x + c) + d
                        if(isempty(x0)), theta0 = optifit.genTheta0(model,xdata,ydata,wts); else theta0 = x0; end
                        fun = @(theta,xdata) theta(1)*sin(theta(2)*xdata + theta(3)) + theta(4);
                        grad = @(theta,xdata) [sin(theta(3) + theta(2)*xdata), theta(1)*xdata.*cos(theta(3) + theta(2)*xdata), theta(1)*cos(theta(3) + theta(2)*xdata), ones(size(xdata))];                            
                        str = 't1*sin(t2*xd + t3) + t4';
                        
                    case 'sin2' %a*sin(b*x + c) + d*sin(e*x + f) + g
                        if(isempty(x0)), theta0 = optifit.genTheta0(model,xdata,ydata,wts); else theta0 = x0; end
                        fun = @(theta,xdata) theta(1)*sin(theta(2)*xdata + theta(3)) + theta(4)*sin(theta(5)*xdata + theta(6)) + theta(7);
                        grad = @(theta,xdata) [sin(theta(3) + theta(2)*xdata), theta(1)*xdata.*cos(theta(3) + theta(2)*xdata), theta(1)*cos(theta(3) + theta(2)*xdata),...
                                               sin(theta(6) + theta(5)*xdata), theta(4)*xdata.*cos(theta(6) + theta(5)*xdata), theta(4)*cos(theta(6) + theta(5)*xdata), ones(size(xdata))];                            
                        str = 't1*sin(t2*xd + t3) + t4*sin(t5*xd + t6) + t7';
                        
                    case 'sin3' %a*sin(b*x + c) + d*sin(e*x + f) + g*sin(h*x + i) + j
                        if(isempty(x0)), theta0 = optifit.genTheta0(model,xdata,ydata,wts); else theta0 = x0; end
                        fun = @(theta,xdata) theta(1)*sin(theta(2)*xdata + theta(3)) + theta(4)*sin(theta(5)*xdata + theta(6)) + theta(7)*sin(theta(8)*xdata + theta(9)) + theta(10);
                        grad = @(theta,xdata) [sin(theta(3) + theta(2)*xdata), theta(1)*xdata.*cos(theta(3) + theta(2)*xdata), theta(1)*cos(theta(3) + theta(2)*xdata),...
                                               sin(theta(6) + theta(5)*xdata), theta(4)*xdata.*cos(theta(6) + theta(5)*xdata), theta(4)*cos(theta(6) + theta(5)*xdata),...
                                               sin(theta(9) + theta(8)*xdata), theta(7)*xdata.*cos(theta(9) + theta(8)*xdata), theta(7)*cos(theta(9) + theta(8)*xdata), ones(size(xdata))];                            
                        str = 't1*sin(t2*xd + t3) + t4*sin(t5*xd + t6) + t7*sin(t8*xd + t9) + t10';    
                        
                    case{'poly11','poly12','poly13','poly21','poly22','poly23','poly31','poly32','poly33'}
                        error('%s is for surface fitting only.\n\nIf you intended to fit a surface, please supply the surface data via optifit(xdata,{ydata,zdata},''%s'').',model,model)
                        
                    otherwise
                        warning(s1);
                        warning(s);
                        error('Unknown model type ''%s''\n',model);                                            
                end  
                warning(s1);
                warning(s);
                %Construct OPTI object and solve
                Opt = opti('fun',fun,'grad',grad,'data',xdata,ydata,'bounds',lb,ub,'x0',theta0,'weights',wts,'options',opts,'name',['OPTIFIT:' str ' [' model ']']);        
                [theta,~,ef] = solve(Opt);
                if(ef ~= 1 && optiWarnLevel(opts.warnings))
                    optiwarn('optifit:nosol','Optimizer did not find a solution. Try a different solver, use a better initial guess, or increase maximum fevals');
                end
                
            %Custom model
            elseif(isa(model,'function_handle'))
                if(isempty(x0))
                    error('You must supply x0 for custom models');
                end
                %Attempt to derive gradient
                if(~isempty(which('syms')))
                    try
                        xs = sym('x',[length(x0),1]);
                        if(nargin(model)==1)
                            jac = jacobian(model(xs),xs);
                            grad = SymBuilder.sym2fun(jac,xs);
                        elseif(nargin(model)==2)
                            jac = jacobian(model(xs,xdata),xs);
                            grad = SymBuilder.sym2fun(jac,xs);
                        else
                            error('Model function should be called as model(x) or model(x,xdata) only!');
                        end
                    catch
                        jac = [];
                        grad = [];
                    end
                else
                    jac = [];
                    grad = [];
                end 
                %If linear in the parameters, guess new x0
                if(~isempty(jac) && isempty(symvar(jac)))
                    x0 = double(jac)\zdata;
                end
                fun = model;
                Opt = opti('fun',model,'grad',grad,'data',xdata,ydata,'bounds',lb,ub,'weights',wts,'x0',x0,'options',opts);
                [theta,~,ef,info] = solve(Opt);
                stats = fitStats(Opt);
                %Try find a better solution with a different solver
                if(ef ~= 1 || stats.Rsquare < 0.8)
                    if(isempty(strfind(info.Algorithm,'NL2SOL')))
                        Opt = opti(Opt,'x0',x0,'solver','NL2SOL');
                        [theta,~,ef] = solve(Opt);
                        stats = fitStats(Opt);                        
                    end
                    if((ef ~= 1 || stats.Rsquare < 0.8) && isempty(strfind(info.Algorithm,'MKLTRNLS')))
                        Opt = opti(Opt,'x0',x0,'solver','MKLTRNLS');
                        [theta,~,ef] = solve(Opt);
                        if(strcmpi(Opt.info.Status,'Exceeded Iterations')) %try with a few more
                            Opt = opti(Opt,'x0',x0,'options',optiset(Opt.opts,'maxiter',3e3));
                            [theta,~,ef] = solve(Opt);
                        end
                        %Try lmder as last chance
                        if(isempty(Opt.prob.lb) && ef~=1)
                            Opt = opti(Opt,'x0',x0,'solver','lmder');
                            [theta,~,ef] = solve(Opt);
                        end
                        try
                            stats = fitStats(Opt);
                        catch
                            optiwarn('optifit:nostats','Could not solve for fit statistics');
                            stats = [];
                        end
                    end
                end
                if(ef ~= 1 && optiWarnLevel(opts.warnings))
                    optiwarn('optifit:nosol','Nonlinear optimizer did not find a solution. Try a different solver, use a better initial guess, or increase maximum fevals');
                end
                
            %Symbolc model
            elseif(isa(model,'sym'))
                error('not implemented yet');
                
            else
                error('Unknown model type ''%s''\n',char(model));
            end         
            
            %Fit Statistics
            if(isempty(stats))
                stats = fitStats(Opt);
            end
        end
        
        %-- Fit 3D Model --%
        function [theta,Opt,stats,fun,grad] = model3DFit(xdata,ydata,zdata,model,x0,lb,ub,wts,opts)
            
            %Standard Models - Use function knowledge to generate reasonable
            %initial guess.
            if(ischar(model))
                switch(lower(model))                    
                    case 'poly11'
                        theta0 = [xdata ydata ones(size(xdata))]\zdata; 
                        fun = @(theta,xdata,ydata) theta(1)*xdata + theta(2)*ydata + theta(3);
                        grad = @(theta,xdata,ydata) [xdata ydata ones(size(xdata))];
                        str = 't1*xd + t2*yd + t3';
                        
                    case 'poly12'
                        theta0 = [ydata.^2 xdata.*ydata xdata ydata ones(size(xdata))]\zdata;
                        fun = @(theta,xdata,ydata) theta(1)*ydata.^2 + theta(2)*xdata.*ydata + theta(3)*xdata + theta(4)*ydata + theta(5);
                        grad = @(theta,xdata,ydata) [ydata.^2 xdata.*ydata xdata ydata ones(size(xdata))];
                        str = 't1*yd^2 + t2*xd*yd + t3*xd + t4*yd + t5';
                        
                    case 'poly13'
                        theta0 = [ydata.^3 ydata.^2 xdata.*ydata.^2 xdata.*ydata xdata ydata ones(size(xdata))]\zdata;
                        fun = @(theta,xdata,ydata) theta(1)*ydata.^3 + theta(2)*ydata.^2 + theta(3)*xdata.*ydata.^2 + theta(4)*xdata.*ydata + theta(5)*xdata + theta(6)*ydata + theta(7);
                        grad = @(theta,xdata,ydata) [ydata.^3 ydata.^2 xdata.*ydata.^2 xdata.*ydata xdata ydata ones(size(xdata))];
                        str = 't1*yd^3 + t2*yd^2 + t3*xd*yd^2 + t4*xd*yd + t5*xd + t6*yd + t7';
                        
                    case 'poly21'
                        theta0 = [xdata.^2 xdata.*ydata xdata ydata ones(size(xdata))]\zdata;
                        fun = @(theta,xdata,ydata) theta(1)*xdata.^2 + theta(2)*xdata.*ydata + theta(3)*xdata + theta(4)*ydata + theta(5);
                        grad = @(theta,xdata,ydata) [xdata.^2 xdata.*ydata xdata ydata ones(size(xdata))];
                        str = 't1*xd^2 + t2*xd*yd + t3*xd + t4*yd + t5';
                        
                    case 'poly22'
                        theta0 = [xdata.^2 ydata.^2 xdata.*ydata xdata ydata ones(size(xdata))]\zdata;
                        fun = @(theta,xdata,ydata) theta(1)*xdata.^2 + theta(2)*ydata.^2 + theta(3)*xdata.*ydata + theta(4)*xdata + theta(5)*ydata + theta(6);
                        grad = @(theta,xdata,ydata) [xdata.^2 ydata.^2 xdata.*ydata xdata ydata ones(size(xdata))];
                        str = 't1*xd^2 + t2*yd^2 + t3*xd*yd + t4*xd + t5*yd + t6';
                        
                    case 'poly23'
                        theta0 = [ydata.^3 xdata.^2 ydata.^2 xdata.^2.*ydata xdata.*ydata.^2 xdata.*ydata xdata ydata ones(size(xdata))]\zdata;
                        fun = @(theta,xdata,ydata) theta(1)*ydata.^3 + theta(2)*xdata.^2 + theta(3)*ydata.^2 + theta(4)*xdata.^2.*ydata + theta(5)*xdata.*ydata.^2 + theta(6)*xdata.*ydata + theta(7)*xdata + theta(8)*ydata + theta(9);
                        grad = @(theta,xdata,ydata) [ydata.^3 xdata.^2 ydata.^2 xdata.^2.*ydata xdata.*ydata.^2 xdata.*ydata xdata ydata ones(size(xdata))];
                        str = 't1*yd^3 + t2*xd^2 + t3*yd^2 + t4*xd^2*yd + t5*xd*yd^2 + t6*xd*yd + t7*xd + t8*yd + t9';
                        
                    case 'poly31'
                        theta0 = [xdata.^3 xdata.^2 xdata.^2.*ydata xdata.*ydata xdata ydata ones(size(xdata))]\zdata;
                        fun = @(theta,xdata,ydata) theta(1)*xdata.^3 + theta(2)*xdata.^2 + theta(3)*xdata.^2.*ydata + theta(4)*xdata.*ydata + theta(5)*xdata + theta(6)*ydata + theta(7);
                        grad = @(theta,xdata,ydata) [xdata.^3 xdata.^2 xdata.^2.*ydata xdata.*ydata xdata ydata ones(size(xdata))];
                        str = 't1*xd^3 + t2*xd^2 + t3*xd^2*yd + t4*xd*yd + t5*xd + t6*yd + t7';
                    
                    case 'poly32'
                        theta0 = [xdata.^3 xdata.^2 ydata.^2 xdata.^2.*ydata xdata.*ydata.^2 xdata.*ydata xdata ydata ones(size(xdata))]\zdata;
                        fun = @(theta,xdata,ydata) theta(1)*xdata.^3 + theta(2)*xdata.^2 + theta(3)*ydata.^2 + theta(4)*xdata.^2.*ydata + theta(5)*xdata.*ydata.^2 + theta(6)*xdata.*ydata + theta(7)*xdata + theta(8)*ydata + theta(9);
                        grad = @(theta,xdata,ydata) [xdata.^3 xdata.^2 ydata.^2 xdata.^2.*ydata xdata.*ydata.^2 xdata.*ydata xdata ydata ones(size(xdata))];
                        str = 't1*xd^3 + t2*xd^2 + t3*yd^2 + t4*xd^2*yd + t5*xd*yd^2 + t6*xd*yd + t7*xd + t8*yd + t9';
                        
                    case 'poly33'
                        theta0 = [xdata.^3 ydata.^3 xdata.^2 ydata.^2 xdata.^2.*ydata xdata.*ydata.^2 xdata.*ydata xdata ydata ones(size(xdata))]\zdata;
                        fun = @(theta,xdata,ydata) theta(1)*xdata.^3 + theta(2)*ydata.^3 + theta(3)*xdata.^2 + theta(4)*ydata.^2 + theta(5)*xdata.^2.*ydata + theta(6)*xdata.*ydata.^2 + theta(7)*xdata.*ydata + theta(8)*xdata + theta(9)*ydata + theta(10);
                        grad = @(theta,xdata,ydata) [xdata.^3 ydata.^3 xdata.^2 ydata.^2 xdata.^2.*ydata xdata.*ydata.^2 xdata.*ydata xdata ydata ones(size(xdata))];
                        str = 't1*xd^3 + t2*yd^3 + t3*xd^2 + t4*yd^2 + t5*xd^2*yd + t6*xd*yd^2 + t7*xd*yd + t8*xd + t9*yd + t10';
                        
                     
                    case{'poly1','poly2','poly3','poly4','poly5','power1','power2','exp1','exp2','exp3','rat01','rat02','rat03','rat11','rat12','rat13','rat21','rat22','rat23',...
                         'sin1','sin2','sin3'}
                        error('%s is for curve fitting only.\n\nIf you intended to fit a curve, please supply the curve data via optifit(xdata,ydata,''%s'').',model,model)
                          
                        
                    otherwise
                        error('Unknown model type ''%s''\n',model);

                end
                ofun = @(theta) fun(theta,xdata,ydata);
                ograd = @(theta) grad(theta,xdata,ydata);
                %Construct OPTI object and solve
                Opt = opti('fun',ofun,'grad',ograd,'ydata',zdata,'bounds',lb,ub,'x0',theta0,'weights',wts,'options',opts,'name',['OPTIFIT:' str ' [' model ']']);        
                [theta,~,ef] = solve(Opt);
                if(ef ~= 1 && optiWarnLevel(opts.warnings))
                    optiwarn('optifit:nosol','Optimizer did not find a solution. Try a different solver, use a better initial guess, or increase maximum fevals');
                end 
                
            %Custom model
            elseif(isa(model,'function_handle')) 
                if(isempty(x0))
                    error('You must supply x0 for custom models');
                end
                %Attempt to derive gradient
                if(~isempty(which('syms')))
                    try
                        xs = sym('x',[length(x0),1]);
                        if(nargin(model)==1)
                            jac = jacobian(model(xs),xs);
                            grad = SymBuilder.sym2fun(jac,xs);
                        elseif(nargin(model)==2)
                            jac = jacobian(model(xs,xdata),xs);
                            grad = SymBuilder.sym2fun(jac,xs);
                        elseif(nargin(model)==3)
                            jac = jacobian(model(xs,xdata,ydata),xs);
                            grad = SymBuilder.sym2fun(jac,xs);
                        else
                            error('Model function should be called as model(x), model(x,xdata) or model(x,xdata,ydata) only!');
                        end
                    catch
                        jac = [];
                        grad = [];
                    end
                else
                    jac = [];
                    grad = [];
                end 
                %Ensure max 1 arg for fcn
                switch(nargin(model))
                    case 2                
                        model = @(x) model(x,xdata);
                    case 3
                        model = @(x) model(x,xdata,ydata);
                end
                %If linear in the parameters, guess new x0
                if(~isempty(jac) && isempty(symvar(jac)))
                    x0 = double(jac)\zdata;
                end
                fun = model;
                Opt = opti('fun',model,'grad',grad,'ydata',zdata,'bounds',lb,ub,'weights',wts,'x0',x0,'options',opts);
                [theta,~,ef,info] = solve(Opt);
                if(ef ~= 1)
                    if(isempty(strfind(info.Algorithm,'NL2SOL')))
                        Opt = opti(Opt,'x0',theta,'solver','NL2SOL');
                        [theta,~,ef] = solve(Opt);
                    elseif(isempty(strfind(info.Algorithm,'MKLTRNLS')))
                        Opt = opti(Opt,'x0',theta,'solver','MKLTRNLS');
                        [theta,~,ef] = solve(Opt);
                    end
                end
                if(ef ~= 1 && optiWarnLevel(opts.warnings))
                    optiwarn('optifit','Nonlinear optimizer did not find a solution. Try a different solver, use a better initial guess, or increase maximum fevals');
                end
                
            end
            
            %Fit Statistics
            stats = fitStats(Opt);                       
        end
           
        
        %-- Generate Initial Guess for Model Types --%
        function theta0 = genTheta0(model,xdata,ydata,wts)
            %TO DO - TAKE WEIGHTS INTO ACCOUNT!
            s = warning('off','MATLAB:polyfit:PolyNotUnique');
            switch(lower(model))
                
                case 'power1'
                    if(any(xdata <= 0)), error('This model requires positive xdata'); end
                    %Use linearization (log10(y) = log10(a) + blog10(x)) == (y = a*x^b))
                    idx = ydata > 0;
                    if(~any(idx)), error('No positive ydata values, try model ''power2'''); end
                    c = polyfit(log10(xdata(idx)),log10(ydata(idx)),1);
                    theta0(1) = 10^(c(2));
                    theta0(2) = c(1);
                    
                case 'power2'
                    if(any(xdata <= 0)), error('This model requires positive xdata'); end
                    %Remove bias then use linearization as above
                    %first check if + or - power
                    tydata = ydata - min(ydata) + sqrt(eps); %ensure positive vals
                    c = polyfit(log10(xdata),log10(tydata),1);
                    if(c(1) < 0) %neg grad
                        bias = ydata(end)-sqrt(eps); %must be better way...
                        idx = ydata-bias > 0;
                        %if we have enough points, skip last few (sensitive to error in bias)                        
                        len = length(ydata(idx));
                        if(len > 3)
                            idx(end-2:end) = false;
                        end
                    else %pos grad
                        bias = ydata(1)-sqrt(eps);
                        idx = ydata-bias > 0;
                        %if we have enough points, skip first few (sensitive to error in bias)                        
                        len = length(ydata(idx));
                        if(len > 3)
                            idx(1:3) = false;
                        end
                    end
                    %Now fit                    
                    c = polyfit(log10(xdata(idx)),log10(ydata(idx)-bias),1);
                    theta0(1) = 10^(c(2));
                    theta0(2) = c(1);
                    theta0(3) = bias;
                
                case 'exp1'
                    %Use linearization (log(y) = log(a) + bx) == (y = a*exp(bx))
                    idx = ydata > 0;
                    if(~any(idx)), error('No positive ydata values, try model ''exp2'''); end
                    if(sum(idx) < 2), error('Not enough positive values to estimate model paramters, try model ''exp2'''); end
                    c = polyfit(xdata(idx),log(ydata(idx)),1);
                    theta0(1) = exp(c(2));
                    theta0(2) = c(1);
                    
                case {'exp2','exp3'}
                    %Remove bias then use linearization as above
                    %first check if + or - exp
                    tydata = ydata - min(ydata) + sqrt(eps); %ensure positive vals
                    c = polyfit(xdata,log(tydata),1);
                    if(c(1) < 0) %neg grad
                        bias = ydata(end)-sqrt(eps); %must be better way...
                        idx = ydata-bias > 0;
                        %if we have enough points, skip last few (sensitive to error in bias)                        
                        len = length(ydata(idx));
                        if(len > 3)
                            idx(end-2:end) = false;
                        end
                    else %pos grad
                        bias = ydata(1)-sqrt(eps);
                        idx = ydata-bias > 0;
                        %if we have enough points, skip first few (sensitive to error in bias)                        
                        len = length(ydata(idx));
                        if(len > 3)
                            idx(1:3) = false;
                        end
                    end
                    %Now fit
                    c = polyfit(xdata(idx),log(ydata(idx)-bias),1);
                    theta0(1) = exp(c(2));
                    theta0(2) = c(1);                    
                    if(strcmpi(model,'exp3'))
                        theta0(3) = exp(c(2))*0.1; %no idea sorry
                        theta0(4) = c(1)*0.1;
                    else
                        theta0(3) = bias;
                    end
                    
                case {'sin1','sin2','sin3'}
                    %Note this assumes uniform spacing between samples!
                    %Also only identifies integer frequencies... See also:
                    %http://www-stat.wharton.upenn.edu/~stine/stat910/lectures/06_harmonic_regr.pdf
                    %ftp://ftp.ti.com/pub/graph-ti/calc-apps/info/sinereg.txt
                    %http://www.math.unl.edu/~sdunbar1/sinefit.html
                    switch(lower(model))
                        case 'sin1', ns = 1;
                        case 'sin2', ns = 2;
                        case 'sin3', ns = 3;
                    end
                    mn = mean(ydata); %find mean
                    %Identify Frequencies
                    fy = fft(ydata-mn); %take fourier transform to identify frequencies
                    [~,idx] = sort(fy(1:floor(length(xdata)/2)),1,'descend'); %only want 1 sided, order for dominant freqs
                    w = 2*pi*(max(0.5,idx(1:ns)-1))./(xdata(end)-xdata(1)); %use index/tspan to find freqs
                    %Construct separable problem (due to double angle trig identity)
                    for i = 1:ns
                        lA(:,i*2-1:i*2) = [sin(w(i).*xdata) cos(w(i).*xdata)];
                    end
                    lf = lA\(ydata-mn); %solve via least squares
                    %Identify Amplitude and Phase
                    theta0=[];
                    for i = 1:ns
                        A = sqrt(lf(i*2-1)^2 + lf(i*2)^2);
                        phi = atan2(lf(i*2),lf(i*2-1));
                        theta0 = [theta0 A w(i) phi]; %#ok<AGROW>
                    end
                    theta0 = [theta0 mn];
                   
            end
            warning(s);
        end                
        
        %-- Attempt to find model structure --%
        function [theta,Opt,stats,fun,grad] = guessModelandFit(xd,yd,x0,lb,ub,wts,opts)
            
            %Model Choosing Parameters
            pgood = 0.1; % preferred maximum pValue for any parameter
            pok = 0.5; % if none satisfy p1, use p2
            rgood = 0.98; % really good fit
            rok = 0.9; % ok fit
            maxA = 3; %max ratio of bounds
            maxWorsePer = 2; % max relative percentage worst than best found still regarded as suitable solution
            
            tops = opts;
            tops.warnings = 'off';
            %Try some quick ones
            d = optifit.modelTestFit([],xd,yd,x0,lb,ub,wts,tops,{'poly1','poly2','poly3','exp1','exp2','power1','power2'});          
            %If nothing is looking (really) good, try a few rationals
            if(~any([d.R] > rgood))
                d = optifit.modelTestFit(d,xd,yd,x0,lb,ub,wts,tops,{'rat01','rat02','rat11','rat12','rat21','rat22'});                
            end
            %If still nothing so far, try some oddball ones
            if(~any([d.R] > rok))
                d = optifit.modelTestFit(d,xd,yd,x0,lb,ub,wts,tops,{'sin1','sin2'});    
            end
            %Concatenate results
            R2 = [d.R]; P = [d.P]; A = [d.A];
            
            %Pursue best AdjR2 so far AND with acceptable max p value           
            idx = P < pgood;
            if(all(~idx))
                idx = P < pok;
                if(all(~idx))
                    idx = true(size(idx)); %give up
                end
            end
            %If more than one suitable, find best based on confidence bound score
            R2 = R2(idx);
            A = A(idx);
            if(length(R2) > 1)
                [r2,sid] = sort(R2,2,'descend');
                a = A(sid);
                if(a(1) > maxA) %assume asymptote or large bounds
                    %find best with A < 3, assuming < 2% worse than best, otherwise accept asymptote/bad bound
                    rOK = abs(r2(2:end) - r2(1))./abs(r2(1))*100 < maxWorsePer & a(2:end) < maxA;
                    if(any(rOK))
                        ind = true(size(R2));
                        i = find(rOK);
                        ind(sid(i(1)+1)) = false; 
                        R2(ind) = 0; %make all the rest zero
                    end
                end
            end            
            [~,ind] = max(R2);
            midx = find(idx); midx = midx(ind);    
            %Index best fit
            theta = d(midx).theta;
            Opt = d(midx).Opt;
            stats = d(midx).s;
            fun = d(midx).f;
            grad = d(midx).g;            
        end
        
        %-- Try a set of models and record results --%
        function d = modelTestFit(d,xd,yd,x0,lb,ub,wts,opts,models)
            
            len = length(d);
            for i = 1:length(models)
                try
                [d(len+i).theta,d(len+i).Opt,d(len+i).s,d(len+i).f,d(len+i).g] = optifit.modelFit(xd,yd,models{i},x0,lb,ub,wts,opts);
                d(len+i).R = d(len+i).s.AdjRsquare;
                d(len+i).P = optifit.pScore(models{i},d(len+i).s.Param.pValues);        
                bnd = d(len+i).s.ConfBnds.bnds(:,2);
                d(len+i).A = abs(max(bnd)./mean(bnd)); %rough check for asymptotes
                catch
                    d(len+i) = struct('theta',[],'Opt',[],'s',[],'f',[],'g',[],'R',0,'P',1,'A',1e6);
                end
            end
        end
        
        %-- Based On Model, Score The pValues --%
        function p = pScore(model,pVals)
            
           switch(lower(model))
               case {'poly1','poly2','poly3','poly4','poly5','power2','exp2','rat01','rat02','rat03'}
                   p = max(pVals(1:end-1)); %ignore bias term
               case {'rat11','rat12','rat13'}
                   idx = true(size(pVals)); idx(2) = false; idx(end) = false;
                   p = max(pVals(idx));
               case {'rat21','rat22','rat23'}
                   idx = true(size(pVals)); idx(3) = false; idx(end) = false;
                   p = max(pVals(idx));
               case 'sin1' %only freq is critical
                   idx = true(size(pVals)); idx([1 3 4]) = false;
                   p = max(pVals(idx));
               case 'sin2'
                   idx = true(size(pVals)); idx([1 3 4 6 7]) = false;
                   p = max(pVals(idx));
               case 'sin3'
                   idx = true(size(pVals)); idx([1 3 4 6 7 9 10]) = false;
                   p = max(pVals(idx));
               otherwise
                   p = max(pVals); %take all
           end
            
        end
    end    
end

