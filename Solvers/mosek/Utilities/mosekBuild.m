function [cmd,param] = mosekBuild(opts)
%MOSEKBUILD  Create MOSEK command string and parameters from options
%
%   [cmd,param] = mosekBuild(opts)

%   This function is based in parts on examples from the MOSEK Toolbox, 
%   Copyright (c) 1998-2011 MOSEK ApS, Denmark.

%   Copyright (C) 2012 Jonathan Currie (I2C2)

%Check we can find mosek
if(exist('mosekopt','file') ~= 3)
    error('It appears MOSEK is not installed on your system');
end

%Get print level information
[plevel,param] = dispLevel(opts.display);

%Get warning status
warn = strcmpi(opts.warnings,'on');

%Create command string
cmd = sprintf('minimize info echo(%d) statuskeys(1) symbcon',plevel);

%MOSEK parameters of MATLAB defined MOSEK options
mset = {'MSK_IPAR_INTPNT_MAX_ITERATIONS','MSK_DPAR_INTPNT_TOL_PFEAS','MSK_DPAR_INTPNT_TOL_DFEAS',...
        'MSK_IPAR_MIO_MAX_NUM_BRANCHES','MSK_IPAR_MIO_MAX_NUM_RELAXS','MSK_DPAR_MIO_TOL_ABS_GAP','MSK_DPAR_MIO_TOL_REL_GAP',...
        'MSK_DPAR_MIO_TOL_ABS_RELAX_INT','MSK_DPAR_OPTIMIZER_MAX_TIME'};
%Options structure fields
mfld = fields(opts);    
%Read and insert MATLAB MOSEK options
for i = 1:length(mset)
    param.(mset{i}) = opts.(mfld{i});
end

%If user has passed their own set of options via mskoption, check then add to param.
if(~isempty(opts.mskoption))
    mkfld = fields(opts.mskoption);
    for i = 1:length(mkfld)
        if(strncmp('MSK_',mkfld{i},4))
            param.(mkfld{i}) = opts.mskoption.(mkfld{i});
        elseif(warn)
            optiwarn('mosek:mskopt','Parameter "%s" does not appear to be a MOSEK option, skipping',mkfld{i});
        end
    end
end

function  [plevel,param] = dispLevel(lev)
%Return MOSEK compatible display level
switch(lower(lev))
    case'off'
        plevel = 0;
        param = struct();
    case 'iter'
        plevel = 3;
        param.MSK_IPAR_LOG_INTPNT = 1;
        param.MSK_IPAR_LOG_SIM = 1;
    case 'final'
        plevel = 3;
        param = struct();
end