function havSym = optiCheckSymTBX()
if(~license('test', 'symbolic_toolbox'))
    optiwarn('OPTI:NoSym','This Function Requires the MATLAB Symbolic Toolbox');
    havSym = false;
else
    havSym = true;
end