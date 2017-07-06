function opti_solverMex(name,src,inc,libs,opts)
%OPTI_SOLVERMEX  Compiles a MEX Interface to an OPTI Solver
%
%   opti_solverMex(name,src,inc,libs,opts)
%
%       src:    Source Files to Compile
%       inc:    Include Directories to Add
%       libs:   Static Libraries to Compile Against
%       opts:   Options Structure with fields below (all optional):
%
%           verb:       True for compiler verbosity (-v)
%           debug:      True for debug build (-g)
%           pp:         Cell array of preprocessors to add (-D)
%           blas:       BLAS/LAPACK libraries to link against ({[]}, MKL, NETLIB)
%           pardiso:    Pardiso Library to link against ({[]}, MKL, Basel)
%           ma57:       MA57 Library to link against ({[]}, MATLAB, HSL)
%           ma27:       MA27 Library to link against ({[]}, HSL)
%           mumps:      Link MUMPS (true/{false}) [v4.0]
%           mumps5:     Link MUMPS v5.0 (true/{false}) [v5.0]
%           asl:        Link AMPL Solver Library (true/{false})
%           expre:      Extra arguments for mex before source file (one string)
%           ifort:      Link against Intel Fortran Dynamic Libraries {false}
%           ifortStatic: Link against Intel Fortran Static Libraries {false}
%           util:       Utility, not a solver, includes extra paths
%           quiet:      Don't print anything

%Process Options
if(nargin > 4)
    if(~isfield(opts,'verb') || isempty(opts.verb)), opts.verb = false; end
    if(~isfield(opts,'debug') || isempty(opts.debug)), opts.debug = false; end
    if(~isfield(opts,'blas')), opts.blas = []; end
    if(~isfield(opts,'pardiso')), opts.pardiso = []; end
    if(~isfield(opts,'ma57')), opts.ma57 = []; end
    if(~isfield(opts,'ma27')), opts.ma27 = []; end
    if(~isfield(opts,'mumps') || isempty(opts.mumps)), opts.mumps = false; end
    if(~isfield(opts,'mumps5') || isempty(opts.mumps5)), opts.mumps5 = false; end
    if(~isfield(opts,'asl') || isempty(opts.asl)), opts.asl = false; end
    if(~isfield(opts,'pp')), opts.pp = []; end
    if(~isfield(opts,'expre')), opts.expre = []; end
    if(~isfield(opts,'util') || isempty(opts.util)), opts.util = false; end
    if(~isfield(opts,'ifort') || isempty(opts.ifort)), opts.ifort = false; end
    if(~isfield(opts,'ifortStatic') || isempty(opts.ifortStatic)), opts.ifortStatic = false; end
    if(~isfield(opts,'quiet') || isempty(opts.quiet)), opts.quiet = false; end
else
    opts.verb = false;
    opts.debug = false;
    opts.blas = [];
    opts.pardiso = [];
    opts.ma57 = [];
    opts.ma27 = [];
    opts.mumps = false;
    opts.mumps5 = false;
    opts.asl = false;
    opts.pp = [];
    opts.expre = [];
    opts.util = false;
    opts.ifort = false;
    opts.ifortStatic = false;
    opts.quiet = false;
end

%Check conflicting BLAS/PARDISO
if(~strcmpi(opts.blas,'mkl') && strcmpi(opts.pardiso,'mkl'))
    if(isempty(opts.blas))
        opts.blas = 'mkl'; %must be used
    else
        error('If using MKL PARDISO you must link against MKL BLAS');
    end
end

%Remove existing mex file from memory
clear(name);

if (~opts.quiet)
    fprintf('\n------------------------------------------------\n');
    fprintf('%s MEX FILE INSTALL\n\n',upper(name));
end

%Build Source File String
if(iscell(src))
    src_str = src{1};
    for i = 2:length(src)
        src_str = sprintf('%s %s',src_str,src{i});
    end
else
    src_str = src;
end

%Build Include String
if(iscell(inc))
    inc_str = [' -I' inc{1}];
    for i = 2:length(inc)
        inc_str = sprintf('%s -I%s',inc_str,inc{i});
    end
elseif(~isempty(inc))
    inc_str = [' -I' inc];
else
    inc_str = '';
end
if(opts.util)
    inc_str = [inc_str ' -I..\..\Solvers\Source\opti '];
else
    inc_str = [inc_str ' -Iopti '];
end

%Build Library String
lib_str = [' -L' getLibPath];
if(opts.util)
    lib_str = [lib_str ' -L..\..\Utilities\Source\' getLibPath];
    lib_str = [lib_str ' -L..\..\Solvers\Source\' getLibPath];
end
if(iscell(libs))
    for i = 1:length(libs)
        lib_str = sprintf('%s -l%s',lib_str,libs{i});
    end
elseif(~isempty(libs))
    lib_str = [lib_str ' -l' libs];
else
    lib_str = '';
end

lib_str = [lib_str ' -llibut '];

%Post Messages (output name + preprocessors)
post = [' -output ' name];
if(isfield(opts,'pp') && ~isempty(opts.pp))
    if(iscell(opts.pp))
        for i = 1:length(opts.pp)
            post = sprintf('%s -D%s',post,opts.pp{i});
        end
    else
        post = [post ' -D' opts.pp];
    end
end

%If compiling with VS2015 but pre R2015b, need to manually add in UCRT location
cc = mex.getCompilerConfigurations();
for i = 1:length(cc)
    if(~isempty(strfind(cc(i).Name,'Microsoft Visual C++')) && str2double(cc(i).Version) >= 14 && verLessThan('matlab','8.6'))
        post = [post opti_FindUCRT()];
        break;
    end
end

%Other Options
if(opts.verb)
    verb = ' -v ';
else
    verb = ' ';
end
if(opts.debug)
    debug = ' -g ';
    post = [post ' -DDEBUG']; %some mex files use this to provide extra info
else
    debug = ' ';
end
if(~isempty(opts.expre))
    opts.expre = [' ' opts.expre ' '];
end

%BLAS Linking
if(isfield(opts,'blas') && ~isempty(opts.blas))
    switch(lower(opts.blas))
        case {'mkl','intel','intelmkl'}
            post = [post ' -DLINK_MKL ' opti_FindMKL];
        case {'mkl_seq','mklseq'}
            post = [post ' -DLINK_MKL ' opti_FindMKL('seq')];
        case 'netlib'
            post = [post ' -DLINK_NETLIB_BLAS -lblas -llapack '];
            opts.ifort = true; %assume compiled with OPTI + Ifort
        otherwise
            error('%s not yet supported for BLAS',opts.blas);
    end
end
%PARDISO Linking
if(isfield(opts,'pardiso') && ~isempty(opts.pardiso))
    switch(lower(opts.pardiso))
        case {'mkl','intel','intelmkl'}
            post = [post ' -DLINK_MKL_PARDISO ']; %mkl linked above as blas
        case {'basel'}
            post = [post ' -DLINK_PARDISO -llibpardiso '];
        otherwise
            error('%s not yet supported for PARDISO',opts.pardiso);
    end
end
%MA57 Linking
if(isfield(opts,'ma57') && ~isempty(opts.ma57))
    switch(lower(opts.ma57))
        case {'hsl'}
            post = [post ' -DLINK_MA57 -llibma57 '];
            post = [post ' -DLINK_METIS -llibmetis ']; %assumed included in build of MA57
            opts.ifort = true; %assume compiled with OPTI + Ifort            
        case {'matlab','ml'}
            post = [post ' -DLINK_ML_MA57 -llibmwma57 '];
        otherwise
            error('%s not yet supported for MA57',opts.ma57);
    end
end
%MA27 Linking
if(isfield(opts,'ma27') && ~isempty(opts.ma27))
    switch(lower(opts.ma27))
        case {'hsl'}
            post = [post ' -DLINK_MA27 -llibma27 '];
        otherwise
            error('%s not yet supported for MA27',opts.ma27);
    end
end
%MUMPS Linking
if(isfield(opts,'mumps') && ~isempty(opts.mumps))
    if(opts.mumps)
        post = [post ' -DLINK_MUMPS -IInclude\Mumps -llibdmumps_c -llibdmumps_f -llibseq_c -llibseq_f -llibmetis -llibpord '];   
        opts.ifort = true;
%         opts.ifortStatic = true; %assume compiled with OPTI + Ifort (statically linked IFort Libs) 
    end
end
%MUMPS v5 Linking
if(isfield(opts,'mumps5') && ~isempty(opts.mumps5))
    if(opts.mumps5)
        post = [post ' -DLINK_MUMPS -IInclude\Mumps5 -llibdmumps5_c -llibdmumps5_f -llibmumps_common_c -llibseq5_c -llibseq5_f -llibpord5 -llibmetis5 '];
        opts.ifort = true; %assume compiled with OPTI + Ifort  
    end
end
%ASL Linking
if(isfield(opts,'asl') && ~isempty(opts.asl))
    if(opts.asl)
        post = [post ' -DLINK_ASL -DNO_STDIO1 -I..\..\Utilities\Source\Include\Asl -L..\..\Utilities\Source\' getLibPath ' -llibasl'];
        if(strcmpi(name,'scip'))
            src_str = [src_str ' scip/ASL/reader_nl.c ']; %include Stefan's ASL reader
        end
    end
end
%Intel Fortran Linking
if(isfield(opts,'ifort') && ~isempty(opts.ifort) && opts.ifort)
    if(isempty(strfind(lower(opts.blas),'mkl'))) %link MKL as well
        [~,fstr,~,~,cmplr] = opti_FindMKL(false,opts.ifortStatic);
        post = [post ' -L"' cmplr '" ' fstr ' '];
    else
        [~,fstr] = opti_FindMKL(false,opts.ifortStatic);
        post = [post ' ' fstr ' '];
    end
end

% Extra preprocessor defines
mver = ver('matlab');
mver = mver.Release;
if (mver(1) == '('), mver = mver(2:end-1); end
post = [post ' -DML_VER=' mver ' -DOPTI_VER='  sprintf('%.2f',optiver)];

%CD to Source Directory
cdir = cd;
if(opts.util)
    cd 'Utilities/Source';
else
    cd 'Solvers/Source';
end
%Compile & Move
try
    evalstr = ['mex' verb debug '-largeArrayDims ' opts.expre src_str inc_str lib_str post];
    if (~opts.quiet)
        fprintf('MEX Call:\n%s\n\n',evalstr);
    end
    eval(evalstr)
    movefile([name '.' mexext],'../','f')
    if (~opts.quiet)
        fprintf('Done!\n');
    end
catch ME
    cd(cdir);
    error('opti:nlopt','Error Compiling %s!\n%s',upper(name),ME.message);
end
cd(cdir);
if (~opts.quiet)
    fprintf('------------------------------------------------\n');
end




function libdir = getLibPath()
% Return architecture dependent general library path
switch(computer)
    case 'PCWIN'
        libdir = 'lib\win32\ ';
    case 'PCWIN64' 
        libdir = 'lib\win64\ ';
    otherwise
        error('This function is only setup for Windows PCs');
end
