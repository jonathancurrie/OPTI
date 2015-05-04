function opti_debug()
% Collects information about your OPTI install, MATLAB and OS and writes to
% a txt file for use in debugging interface problems.

fid = fopen('opti_debug.txt','w');
if(fid < 0)
    error('Cannot create a debug file in the current directory. Please re run in a directory you have write permissions');
end

try
    %Header
    fprintf(fid,'OPTI Debug File generated %s\n\n',datestr(now));

    %Environment Variable Information
    fprintf(fid,'OS: %s\n',getenv('OS'));
    fprintf(fid,'ARCH: %s\n',getenv('PROCESSOR_ARCHITECTURE'));
    fprintf(fid,'IDENT: %s\n',getenv('PROCESSOR_IDENTIFIER'));
    fprintf(fid,'LEVEL: %s\n',getenv('PROCESSOR_LEVEL'));
    fprintf(fid,'REV: %s\n',getenv('PROCESSOR_REVISION'));
    fprintf(fid,'#PROC: %s\n',getenv('NUMBER_OF_PROCESSORS'));

    %MATLAB Information
    mver = ver('MATLAB');
    fprintf(fid,'\nMATLAB: %s, v%s, %s\n',mver.Release,mver.Version,mver.Date);
    fprintf(fid,'MATLAB OS: %s\n',computer);
    
    %Memory Information
    [UV,SV] = memory();
    fprintf(fid,'\nTotal RAM: %d MB\n',int32(SV.PhysicalMemory.Total/1e6));
    fprintf(fid,'Available RAM: %d MB\n',int32(SV.PhysicalMemory.Available/1e6));
    fprintf(fid,'MATLAB Memory: %d MB\n',int32(UV.MemUsedMATLAB/1e6));

    %Pre Req Information
    [havVC,havIC,havIF] = opti_PreReqCheck(false);
    havVCs = 'YES'; havICs = 'YES'; havIFs = 'YES';
    if(~havVC), havVCs = 'NO'; end
    if(~havIC), havICs = 'NO'; end
    if(~havIF), havIFs = 'NO'; end
    fprintf(fid,'\nPreReq Availability:\n');
    fprintf(fid,'Microsoft VC++: %s\n',havVCs);
    fprintf(fid,'Intel C++: %s\n',havICs);
    fprintf(fid,'Intel Fortran: %s\n\n',havIFs);
    
    %MEX Information
    cc = mex.getCompilerConfigurations();
    try
        fprintf(fid,'MEX: %s, ver %s\n',cc.Name,cc.Version);
    catch
        fprintf(fid,'MEX: Not Setup\n');
    end  
    
    %MKL Information
    try
        [~,~,~,mkl_lib,~,mkl_ver] = opti_FindMKL();
        fprintf(fid,'MKL: %s, ver %s\n',mkl_lib,mkl_ver);
    catch
        fprintf(fid,'MKL: Not Found\n');
    end
    
    %CPLEX Information
    try
        [~,~,cplx_lib,~,cplx_ver] = opti_FindCplex();
        fprintf(fid,'CPLEX: %s, ver %s\n',cplx_lib,cplx_ver);
    catch
        fprintf(fid,'CPLEX: Not Found\n');
    end
    
    %Now get OPTI information
    p = which('optiver.m');
    if(~isempty(p))
        %Print Version
        if(optiSolver('scip',0)), str = 'Academic'; else str = 'Open Source'; end
        fprintf(fid,'\nOPTI: v%4.2f (%s)\n',optiver,str);
        %Print Install Location
        fprintf(fid,'Install Location: %s\n',which('opti_Install.m'));
        
        %Print Solver Availability
        fprintf(fid,'\nSolver Availability:\n');
        s = optiSolver('ver');
        for i = 1:length(s)
            fprintf(fid,'%s\n',s{i});
        end
    else
       fprintf('OPTI: Cannot find optiver, OPTI not installed?\n\n'); 
    end
    
    %Tell user all done
    fprintf('Debug information written to ''opti_debug.txt'' in %s\n',cd);
catch ME
    fclose(fid);
    rethrow(ME);
end

fclose(fid);

