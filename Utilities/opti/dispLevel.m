function level = dispLevel(disp_str,finalLev,iterLev)
%DISPLEVEL  Convert String display format to numerical

if(isnumeric(disp_str))
    level = disp_str;
end

if(nargin < 3), iterLev = 2; end
if(nargin < 2), finalLev = 1; end

switch(disp_str)
    case 'off'
        level = 0;
    case 'final'
        level = finalLev;
    case 'iter'
        level = iterLev;
end

end

