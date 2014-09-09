function havSym = optiCheckSymTBX()
c = which('mupadmex');
if(isempty(c))
    optiwarn('OPTI:NoSym','This Function Requires the MATLAB Symbolic Toolbox');
    havSym = false;
else
    havSym = true;
end