function str = optiMLVer()

v = ver('matlab');
if(v.Release(1) == '(')
    str = v.Release(2:end-1);
else
    str = v.Release;
end