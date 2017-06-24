@echo off
rem MSVC141OPTS.BAT
rem
rem    Compile and link options used for building MEX-files
rem    using the Microsoft Visual C++ compiler version 15.0
rem
rem StorageVersion: 1.0
rem C++keyFileName: MSVC141OPTS.BAT
rem C++keyName: Microsoft Visual C++ 2017
rem C++keyManufacturer: Microsoft
rem C++keyVersion: 15.0
rem C++keyLanguage: C++
rem C++keyLinkerName: Microsoft Visual C++ 2017
rem C++keyLinkerVer: 15.0
rem
rem ********************************************************************
rem General parameters
rem ********************************************************************

set MATLAB=%MATLAB%
set VSINSTALLDIR=C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Tools\MSVC\14.10.25017
set VCINSTALLDIR=%VSINSTALLDIR%
rem In this case, LINKERDIR is being used to specify the location of the SDK
set LINKERDIR=C:\Program Files (x86)\Windows Kits\10\bin\10.0.15063.0
set LINKERBIN=%LINKERDIR%
rem Split the linker path into folder and version number
For %%A in ("%LINKERDIR%") do (
    Set LINKERPATH=%%~dpA
    Set LINKERVER="%%~nxA"
)
set LINKERINC=C:\Program Files (x86)\Windows Kits\10\Include\10.0.15063.0
set LINKERLIB=C:\Program Files (x86)\Windows Kits\10\Lib\10.0.15063.0
set PATH=%VCINSTALLDIR%\bin\HostX64\x64;%VSINSTALLDIR%\..\..\..\..\Common7\IDE;%VSINSTALLDIR%\..\..\..\..\Common7\Tools;%LINKERBIN%\x64;%MATLAB_BIN%;%PATH%
set INCLUDE=%VCINSTALLDIR%\INCLUDE;%VCINSTALLDIR%\ATLMFC\INCLUDE;%LINKERINC%\um;%LINKERINC%\shared;%LINKERINC%\winrt;%INCLUDE%
set LIB=%VCINSTALLDIR%\LIB\x64;%VCINSTALLDIR%\ATLMFC\LIB\x64;%LINKERLIB%\um\x64;%MATLAB%\extern\lib\win64;%LIB%
set MW_TARGET_ARCH=win64
echo %LINKERLIBUM%
echo %LIB%

rem ********************************************************************
rem Compiler parameters
rem ********************************************************************
set COMPILER=cl
set COMPFLAGS=/c /GR /W3 /EHs /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_DEPRECATE /D_SECURE_SCL=0 /DMATLAB_MEX_FILE /nologo /MD
set OPTIMFLAGS=/O2 /Oy- /DNDEBUG
set DEBUGFLAGS=/Z7
set NAME_OBJECT=/Fo

rem ********************************************************************
rem Linker parameters
rem ********************************************************************
set LIBLOC=%MATLAB%\extern\lib\win64\microsoft
set LINKER=link
set LINKFLAGS=/dll /export:%ENTRYPOINT% /LIBPATH:"%LIBLOC%" libmx.lib libmex.lib libmat.lib /MACHINE:X64 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /manifest /incremental:NO /implib:"%LIB_NAME%.x" /MAP:"%OUTDIR%%MEX_NAME%%MEX_EXT%.map"
set LINKOPTIMFLAGS=
set LINKDEBUGFLAGS=/debug /PDB:"%OUTDIR%%MEX_NAME%%MEX_EXT%.pdb"
set LINK_FILE=
set LINK_LIB=
set NAME_OUTPUT=/out:"%OUTDIR%%MEX_NAME%%MEX_EXT%"
set RSP_FILE_INDICATOR=@

rem ********************************************************************
rem Resource compiler parameters
rem ********************************************************************
set RC_COMPILER=rc /fo "%OUTDIR%mexversion.res"
set RC_LINKER=

set POSTLINK_CMDS=del "%LIB_NAME%.x" "%LIB_NAME%.exp"
set POSTLINK_CMDS1=mt -outputresource:"%OUTDIR%%MEX_NAME%%MEX_EXT%;2" -manifest "%OUTDIR%%MEX_NAME%%MEX_EXT%.manifest"
set POSTLINK_CMDS2=del "%OUTDIR%%MEX_NAME%%MEX_EXT%.manifest"
set POSTLINK_CMDS3=del "%OUTDIR%%MEX_NAME%%MEX_EXT%.map"
