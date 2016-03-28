function [ucrt_link,ucrt_inc,ucrt_lib,ucrt_ver] = opti_FindUCRT()
%Finds Windows Universal C Runtime installation folder

%Known UCRT path locations (Modify to suit your system by adding to cell arrays, or create a new structure for other versions)
UCRT10 = {'C:\Program Files (x86)\Windows Kits\10'};

%Add any new structures to below
UCRTLIB = {UCRT10};


%DO NOT MODIFY BELOW HERE
for i = 1:length(UCRTLIB)
    %Check for MKL
    [ok,ucrt_link,ucrt_inc,ucrt_lib,ucrt_ver] = checkUCRTVer(UCRTLIB{i}); 
    if(ok), return; end
end

%Still no Luck - user must locate
error('Could not find the Universal CRT Location on your computer. Please modify this file to locate it');


function [ok,ucrt_link,ucrt_inc,ucrt_lib,ucrt_ver] = checkUCRTVer(ucrtstr)
%Local check function
ok = false; ucrt_link ='';
%Find MKL
for i = 1:length(ucrtstr)
    if(exist(ucrtstr{i},'dir'))
        dir_ucrt = ucrtstr{i};
        %Check for Include Directory
        ucrt_inc = [dir_ucrt '\Include'];
        if(~exist(ucrt_inc,'dir'))
            return; %no luck
        end
        %Find highest version of UCRT
        vers = dir(ucrt_inc);
        maxVer = 0; idx = 0;
        for j = 1:length(vers)
            if(length(vers(j).name) > 3)
                try
                    tVer = str2double(strrep(vers(j).name,'.',''));
                    if(tVer > maxVer)
                        maxVer = tVer;
                        idx = j;
                    end
                catch
                end
            end
        end
        if(idx > 0)
            ucrt_ver = vers(idx).name;
        else
            return; %no luck
        end
        %Assign include directory
        ucrt_inc = [ucrt_inc filesep ucrt_ver filesep 'ucrt']; 
        if(~exist(ucrt_inc,'dir'))
            return;
        end            
        %Get Library Directory
        switch(computer)
            case 'PCWIN'
                ucrt_lib = [dir_ucrt filesep 'Lib' filesep ucrt_ver filesep 'ucrt\x86'];
                ucrt_link = ' -lucrt';
            case 'PCWIN64'
                ucrt_lib = [dir_ucrt filesep 'Lib' filesep ucrt_ver filesep 'ucrt\x64'];  
                ucrt_link = ' -lucrt';
        end
        if(~exist(ucrt_lib,'dir'))
            return;
        end    
        %Build complete linker string
        ucrt_link = [' -I"' ucrt_inc '" -L"' ucrt_lib '" ' ucrt_link ];   
        ok = true;
        return;
    end
end
