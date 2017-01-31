function stop = optiplotfval(iter,fval,~)
%OPTIPLOTVAL  Plot objective vs iterations
%
%   Default plotting callback for scalar objectives
%
%   Based on MATLAB's optimplotfval

persistent liter

%Ensure always create a new plot
if(isempty(liter))
    liter = Inf;
end

%Never stop from this function
stop = false;

%Only plot scalars
if(length(fval) > 1)
    fval = norm(fval);
end

%Check for new solver run
if(iter < liter)
    h = plot(iter,fval,'b*');
    title(sprintf('Current Function Value: %g',fval));
    xlabel('Iteration /  Function Evaluation');
    set(h,'Tag','optiplotf1');
    ylabel('Function Value')
else
    h = findobj(get(gca,'Children'),'Tag','optiplotf1');
    newX = [get(h,'Xdata') iter];
    newY = [get(h,'Ydata') fval];
    set(h,'Xdata',newX, 'Ydata',newY);
    set(get(gca,'Title'),'String',sprintf('Current Function Value: %g',fval));
end

%Save for future comparison
liter = iter;

%Call drawnow to refresh the screen
drawnow;