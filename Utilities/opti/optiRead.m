function prob = optiRead(filename)
%optiRead  Read a Mathematical File and convert to Matlab Values
%
% prob = optiRead(filename) reads the file specified by filename and
% guesses the problem type by the file extension, then converts it to an
% optiprob. The returned structure is solver independent so the user can 
% manually extract matrices as required or supply it directly to opti().
%
%   Available File Types:
%       - MPS  [.mps]
%       - QPS  [.qps]
%       - LP   [.lp]        (Appears to be CPLEX LP version - Simplified)
%       - AMPL [.mod,.nl]   (Note GMPL files are not recognised with this interface)
%       - SDPA-S [.dat-s]   (Sparse Format)
%       - SDPA   [.dat]     (Dense Format)
%       - SeDuMi [.mat]     
%
% You may specify a full path to the file, or if you specify a filename
% only, it must exist on the MATLAB path.
%
% The routines underneath use COIN-OR & GLPK utilities for File IO. See 
% attached EPL License for COIN-OR and GPL for GLPK. In addition, it uses 
% the Netlib AMPL Solver Library code. You must have a licensed version of 
% AMPL (or the student edition) present on your computer to read .mod files. 
% For more information see www.ampl.com.

%   Copyright (C) 2013 Jonathan Currie (Control Engineering)


if(isempty(strfind(filename,'.')))
    error('You must supply a file extension to this function!');
end
ext = regexp(filename,'\.','split');
switch(lower(ext{end}))
    case {'dat','mat','dat-s'}
        prob = sdpRead(filename);
    case {'mps','qps','lp'}            
        prob = coinRead(filename);
    case {'nl','mod'}
        prob = amplRead(filename);
    otherwise
        error('Unknown file extension to read from!');
end
