function [havVC,havIC,havIF] = opti_PreReqCheck(doCD)

if(nargin < 1), doCD = true; end %assumes in base OPTI directory!

havVC = true;
havIC = true;
havIF = true;

%Check for VC++ 2013
if(doCD), cd('../../Solvers/'); end
try
    a = nlopt; %#ok<NASGU>
catch
    havVC = false;
end
%Check for IC 2015 [not req from OPTI v >= 2.12]
try
    a = clp; %#ok<NASGU>
catch
    havIC = false;
end
%Check for IFort 2015
try
    a = lmder; %#ok<NASGU>
catch
    havIF = false;
end


%% OLD CODE
% function str = opti_PreReqCheck(mode,verb)
% %Check for VC++ 2010 / VC++ 2012 on x86 and x64 systems
% 
% if(nargin < 2), verb = 1; end
% if(nargin < 1), mode = 'VS2012'; end
% 
% ROOTKEY = 'HKEY_LOCAL_MACHINE';
% switch(lower(mode))
%     case {'2012','vs2012'}
%         %VC++ 2012 Keys
%         KEY32 = 'SOFTWARE\Wow6432Node\Microsoft\DevDiv\VC\Servicing\11.0\RuntimeMinimum';
%         KEY64 = 'SOFTWARE\Microsoft\DevDiv\VC\Servicing\11.0\RuntimeMinimum1';
%         rtstr = 'VC++ 2012';
%         rtredist86 = '';
%         rtrunt86 = '';
%         rtredist64 = '';
%         rtrunt64 = '';
%         instkey = 'Install';
%     case {'2010','vs2010'}
%         %VC++ 2010 Keys
%         KEY32 = 'SOFTWARE\Microsoft\VisualStudio\10.0\VC\';
%         KEY64 = 'SOFTWARE\Wow6432Node\Microsoft\VisualStudio\10.0\VC\';
%         rtstr = 'VC++ 2010';
%         rtredist86 = 'VCRedist\x86';
%         rtrunt86 = 'Runtimes\x86';
%         rtredist64 = 'VCRedist\x64';
%         rtrunt64 = 'Runtimes\x64';
%         instkey = 'Installed';
%     otherwise
%         error('Unknown Visual C++ / Intel C++ Version to check for! Only VC++ 2010 and 2012 supported');
% end
% 
% mver = ver('MATLAB');
% if(verb), fprintf('- Checking operating system...\n'); end
% switch(mexext)
%     case 'mexw32'
%         if(verb), fprintf('MATLAB %s 32bit (Windows x86) detected\n',mver.Release); end
%     case 'mexw64'
%         if(verb), fprintf('MATLAB %s 64bit (Windows x64) detected\n',mver.Release); end
%     otherwise
%         error('OPTI Toolbox is compiled only for Windows systems - sorry!');
% end
% 
% if(verb), fprintf('\n- Checking for required pre-requisites...\n'); end
% 
% str = [];
% switch(mexext)
%     %If mexw32 - could be 32bit windows or 64bit windows with 32bit matlab, must check both!
%     case 'mexw32'
%         %Check 32bit Location, Redist
%         try
%             if(winqueryreg(ROOTKEY,[KEY32 rtredist86],instkey))
%                 str = [rtstr ' x86 Redistributable (Win32) Found'];
%             end
%         catch %#ok<*CTCH>
%             %Check 32bit Location, Runtime
%             try
%                 if(winqueryreg(ROOTKEY,[KEY32 rtrunt86],instkey))
%                     str = [rtstr ' x86 Runtime (Win32) Found'];
%                 end
%             catch
%                 %Check 64bit Location, Redist
%                 try
%                     if(winqueryreg(ROOTKEY,[KEY64 rtredist86],instkey))
%                         str = [rtstr ' x86 Redistributable (Win64) Found'];
%                     end
%                 catch
%                     %Check 64bit Location, Runtime
%                     try
%                         if(winqueryreg(ROOTKEY,[KEY64 rtrunt86],instkey))
%                             str = [rtstr ' x86 Redistributable (Win64) Found'];
%                         end
%                     catch
%                         error(prereqerror('x86',rtstr));
%                     end
%                 end
%             end
%         end                                    
%     case 'mexw64'
%         %Check Redist
%         try
%             if(winqueryreg(ROOTKEY,[KEY64 rtredist64],instkey))
%                 str = [rtstr ' x64 Redistributable Found'];
%             end
%         catch
%             %Check Runtime
%             try
%                 if(winqueryreg(ROOTKEY,[KEY64 rtrunt64],instkey))
%                     str = [rtstr ' x64 Runtime Found'];
%                 end                
%             catch
%                 error(prereqerror('x64',rtstr));
%             end            
%         end
% end
% 
% if(isempty(str))
%     error('Could not determine the required pre-requisites, ensure you are running Windows x86 or x64!');
% else
%     if(verb), fprintf('%s\n',str); end
% end
% 
% 
% function str = prereqerror(build,rtstr)
% 
% switch(rtstr)
%     case 'VC++ 2012'
%         l = 'http://www.microsoft.com/en-us/download/details.aspx?id=30679';
%         m = 'Microsoft';
%     case 'VC++ 2010'
%         switch(build)
%             case 'x64'
%                 l = 'http://www.microsoft.com/download/en/details.aspx?id=14632';
%             otherwise
%                 l = 'http://www.microsoft.com/download/en/details.aspx?id=5555';
%         end
%         m = 'Microsoft';
%     case 'Intel C++ XE 2013'
%         l = 'http://software.intel.com/en-us/articles/redistributable-libraries-for-intel-c-and-visual-fortran-composer-xe-2013-for-windows';
%         m = 'Intel';
% end
%         
% str = [sprintf(['Cannot find the %s ' rtstr ' ' build ' Redistributable!\n\n'],m)...
% 	   sprintf('Please close MATLAB and download the redistributable from %s:\n',m)...
%        l...
%        sprintf('\n\nThen reinstall OPTI Toolbox.')];