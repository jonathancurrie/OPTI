function optiPrintFitStats(fitstats,x)
%OPTIPRINTFITSTATS  Print OPTI Regression Fit Statistics

disp('---------------------------------------------------------------------------------'); 
fprintf('OPTI Fit Statistics\n');
fprintf(' R^2:             %g\n',fitstats.Rsquare);
fprintf(' Adjusted R^2:    %g\n',fitstats.AdjRsquare);
fprintf(' RMSE:            %g\n',fitstats.RMSE);
fprintf(' SSE:             %g\n',fitstats.SSE);
if(isfield(fitstats,'ModelStructure') && ~isempty(fitstats.ModelStructure))
    fprintf(' Model Structure: %s\n',fitstats.ModelStructure);
    var = 't';
else
    var = 'x';
end
if(fitstats.Model.Int)
    fprintf(' Model Form:      Constant (With Intercept)\n');
else
    fprintf(' Model Form:      Zero (No Intercept)\n');    
end
if(strcmpi(fitstats.SolverStatus,'Exceeded Iterations') || ~isempty(strfind(fitstats.SolverStatus,'Error')) || ~isempty(strfind(fitstats.SolverStatus,'Infeasible')))
    fprintf(2,' Solver Status:   ''%s''\n',fitstats.SolverStatus);
else
    fprintf(' Solver Status:   ''%s''\n',fitstats.SolverStatus);
end

nparam = length(x);

%Print ANOVA Information
fprintf('\nNonlinear Least-Squares Analysis Of Variance:\n');

fprintf(' Source          DF     Sum of Squares      Mean Square       F Value       p Value\n');
fprintf(' Model         %4d       %12g     %12g  %12g  %12g\n',fitstats.DF(1),fitstats.SOS(1),fitstats.MS(1),fitstats.Model.FStat,fitstats.Model.pValue);
fprintf(' Error         %4d       %12g     %12g\n',fitstats.DF(2),fitstats.SOS(2),fitstats.MS(2));
fprintf(' Uncorrected   %4d       %12g\n',fitstats.DF(3),fitstats.SOS(3));

%Print Confidence Information
fprintf('\nNonlinear Least-Squares Confidence Interval & Coefficient Statistics:\n');

fprintf(' Parameter      Estimate (%5g%% CI)     Std Error       t Value       p Value\n',fitstats.Conf*100);
for i = 1:nparam
    try
        if(isnan(fitstats.ConfInt(i)))
            if(fitstats.BndIdx(i))
                fprintf('  %s%-3d  %12g  (%12s)  %12s  %12s  %12s\n',var,i,x(i),'At Bound  ','-','-','-');
            else
                fprintf('  %s%-3d  %12g  (%12s)  %12s  %12s  %12s\n',var,i,x(i),'NaN/Inf  ','-','-','-');
            end
        else
            fprintf('  %s%-3d  %12g ±(%12.6g)  %12g  %12g  %12g\n',var,i,x(i),fitstats.ConfInt(i),fitstats.Param.StdError(i),fitstats.Param.tStat(i),fitstats.Param.pValues(i));
        end
    catch
        fprintf(2,'  %s%-3d          (Error)\n',var,i);
    end
end

disp('---------------------------------------------------------------------------------'); 