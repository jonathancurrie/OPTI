function x0 = movex0(lb,ub,x0,warn)
%MOVEX0 Ensure x0 is within problem bounds

%   Copyright (C) 2011 Jonathan Currie (IPL)

if(~exist('warn','var') || isempty(warn))
    warn = 1;
end

if(~isempty(lb))
    lbind = x0 < lb;
    if(any(lbind))
        x0(lbind) = lb(lbind);
    end
else
    lbind = 0;
end
if(~isempty(ub))
    ubind = x0 > ub;
    if(any(ubind))
        x0(ubind) = ub(ubind);
    end
else
    ubind = 0;
end

if((warn>1) && (any(lbind) || any(ubind)))
    optiwarn('opti:x0','x0 was not within the problem bounds. It has been moved to suit NLP / NLS interface so it is within the bounds.');
end

end

