function libdir = opti_GetLibPath()
% Return architecture dependent general library path

switch(computer)
    case 'PCWIN'
        libdir = 'lib\win32\ ';
    case 'PCWIN64' 
        libdir = 'lib\win64\ ';
    otherwise
        error('This function is only setup for Windows PCs');
end

end

