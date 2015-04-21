function [mkl_link,mkl_forstr,mkl_inc,mkl_lib,mkl_cmplr,mkl_ver] = opti_FindMKL(seq)
%Finds Intel MKL Libraries and returns path and version

%Sequential build stuff
if(nargin && ((ischar(seq) && strcmpi(seq,'seq')) || seq))
    seq = true;
else
    seq = false;
end

%Known MKL path locations (Modify to suit your system by adding to cell arrays, or create a new structure for other versions)
MKL11_2.compiler = {'C:\Program Files\Intel\Composer XE 2015\compiler','C:\Program Files (x86)\Intel\Composer XE 2015\compiler'};
MKL11_2.mkl = {'C:\Program Files\Intel\Composer XE 2015\mkl','C:\Program Files (x86)\Intel\Composer XE 2015\mkl'};
MKL11_2.ver = '11.2';

MKL11_1.compiler = {'C:\Program Files\Intel\Composer XE 2013 SP1\compiler','C:\Program Files (x86)\Intel\Composer XE 2013 SP1\compiler'};
MKL11_1.mkl = {'C:\Program Files\Intel\Composer XE 2013 SP1\mkl','C:\Program Files (x86)\Intel\Composer XE 2013 SP1\mkl'};
MKL11_1.ver = '11.1';

MKL11.compiler = {'C:\Program Files\Intel\Composer XE 2013\compiler','C:\Program Files (x86)\Intel\Composer XE 2013\compiler','D:\Program Files (x86)\Intel\Composer XE 2013\compiler'};
MKL11.mkl = {'C:\Program Files\Intel\Composer XE 2013\mkl','C:\Program Files (x86)\Intel\Composer XE 2013\mkl','D:\Program Files (x86)\Intel\Composer XE 2013\mkl'};
MKL11.ver = '11';

MKL10_3_1.compiler = {'C:\Program Files\Intel\Composer XE 2011 SP1\compiler\','C:\Program Files (x86)\Intel\Composer XE 2011 SP1\compiler\'};
MKL10_3_1.mkl = {'C:\Program Files\Intel\Composer XE 2011 SP1\mkl\','C:\Program Files (x86)\Intel\Composer XE 2011 SP1\mkl\'};
MKL10_3_1.ver = '10.3 SP1';

MKL10_3.compiler = {'C:\Program Files\Intel\ComposerXE-2011\compiler\','C:\Program Files (x86)\Intel\ComposerXE-2011\compiler\'};
MKL10_3.mkl = {'C:\Program Files\Intel\ComposerXE-2011\mkl\','C:\Program Files (x86)\Intel\ComposerXE-2011\mkl\'};
MKL10_3.ver = '10.3';

%Add any new structures to below
MKLLIB = {MKL11_2,MKL11_1,MKL11,MKL10_3_1,MKL10_3};


%DO NOT MODIFY BELOW HERE
for i = 1:length(MKLLIB)
    %Check for MKL
    [ok,mkl_link,mkl_forstr,mkl_inc,mkl_lib,mkl_cmplr,mkl_ver] = checkMKLVer(MKLLIB{i},seq); 
    if(ok), return; end
end

%Still no Luck - user must locate
error('Could not find the Intel MKL Location on your computer. Please modify this file to locate it');


function [ok,mkl_link,mkl_forstr,mkl_inc,mkl_lib,mkl_cmplr,mkl_ver] = checkMKLVer(mklstr,seq)
%Local check function
ok = false; mkl_link ='';
%Find MKL
for i = 1:length(mklstr.mkl)
    if(exist(mklstr.mkl{i},'dir'))
        dir_mkl = mklstr.mkl{i};
        %Get Include Directory
        mkl_inc = [dir_mkl '\include'];
        if(~exist(mkl_inc,'dir'))
            error('Could not find MKL include directory. Checked:\n %s',mkl_inc);
        end
        %Get Library Directory
        switch(computer)
            case 'PCWIN'
                mkl_lib = [dir_mkl '\lib\ia32'];
                mkl_mklstr = ' -lmkl_intel_c';
            case 'PCWIN64'
                mkl_lib = [dir_mkl '\lib\intel64'];  
                mkl_mklstr = ' -lmkl_intel_lp64';
        end
        if(~exist(mkl_lib,'dir'))
            error('Could not find MKL lib directory. Checked:\n %s\n',mkl_lib);
        end    
        %Complete lib string
        if(seq)
            mkl_mklstr = [mkl_mklstr ' -lmkl_sequential -lmkl_core ']; %#ok<AGROW>
        else
            mkl_mklstr = [mkl_mklstr ' -lmkl_intel_thread -lmkl_core ']; %#ok<AGROW>
        end
        %Build complete linker string
        mkl_link = [' -I"' mkl_inc '" -L"' mkl_lib '" ' mkl_mklstr ];
        mkl_ver = mklstr.ver;   
        ok = true;
        break;
    end
end

%Find Compiler
for i = 1:length(mklstr.compiler)
    if(exist(mklstr.compiler{i},'dir'))
        dir_cmplr = mklstr.compiler{i};
        %Get Library Directory
        switch(computer)
            case 'PCWIN'
                mkl_cmplr = [dir_cmplr '\lib\ia32'];
            case 'PCWIN64'
                mkl_cmplr = [dir_cmplr '\lib\intel64'];                  
        end
        if(~exist(mkl_cmplr,'dir'))
            error('Could not find MKL compiler directory. Checked:\n %s\n',mkl_cmplr);
        end  
        %Complete linker string
        mkl_link = [mkl_link ' -L"' mkl_cmplr '" -llibiomp5md ']; %#ok<AGROW>
        %Optional Fortran linker string
        mkl_forstr = ' -lifconsol -llibifcoremd -llibifportmd -llibmmd -llibirc -lsvml_disp -lsvml_dispmd ';
        break;
    end
end

if(isempty(mkl_link))
    error('Could Not Find Intel MKL On Your System! Modify this file to locate it on your PC.');
end
