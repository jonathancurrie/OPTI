function level = optiWarnLevel(warn)
%Return opti warning level in numerical form

%Get Warning Level
if(strcmpi(warn,'all'))
    level = 2;
elseif(strcmpi(warn,'critical'))
    level = 1;
else
    level = 0;
end