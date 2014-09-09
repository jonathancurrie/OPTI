function h = plotConfReg(confStats)
%PLOTCONFREG  Plot Confidence Bounds on Data Fitting Plot

std2c = [250 213 176]/255;
lgray = [0.6 0.6 0.6]*1.3;
nc = size(confStats.ConfBnds.bnds,2)/2;
x = confStats.ConfBnds.xdata;
%If matrix, skip plotting
if(size(x,1) > 1 && size(x,2) > 1)
    optiwarn('opti:confxdatamat','Don''t know how to plot confidence bounds with xdata as a matrix...');
    h = [];
    return;
end
if(all(isnan(x)))
    h = [];
    return;
end
if(size(x,2) > 1), x = x'; end
try
    h(1) = patch([x;flipud(x)],[confStats.ConfBnds.bnds(:,1);flipud(confStats.ConfBnds.bnds(:,2))],std2c,'edgecolor',lgray);

    for i = 2:nc
        patch([x;flipud(x)],[confStats.ConfBnds.bnds(:,i*2-1);flipud(confStats.ConfBnds.bnds(:,i*2))],std2c,'edgecolor',lgray);
    end
catch ME
    optiwarn('opti:confxdatamat','Error plotting confidence region (sorry)...\n %s',ME.message);
    h = [];
    return;
end
    