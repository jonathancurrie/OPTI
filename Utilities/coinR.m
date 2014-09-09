% COINR  Read a MPS, QPS, GMPL, GAMS or LP file into Matlab
%
% THIS IS A LOW LEVEL FUNCTION - USE coinRead() INSTEAD!
%
% coinR uses the COIN Utilities library to open and parse the file. GLPK is
% also used to parse GMPL files. 
%
%   [f,A,rl,ru,lb,ub,ivars,H,name,sostype,sosind,soswt,objc] = coinR(path,filetype,print) 
%
%   Input arguments:
%       path - full path to the file on your PC.
%       filetype - one of the following types: [mps,qps,lp,mod,gms]
%       print - 1 to enable debug printing, 0 to disable
%
%   Return arguments:
%       f - linear objective vector
%       A - linear inequality + equality matrix (sparse)
%       rl - linear constraint lower bounds
%       ru - linear constraint upper bounds
%       lb - decision variable lower bounds
%       ub - decision variable upper bounds
%       ivars - integer variable indicies
%       H - quadratic objective matrix (sparse)
%       name - problem name
%       sostype - SOS types
%       sosind - SOS indices
%       soswt - SOS weights
%       objbias - Objective bias term
%
%
%   Copyright (C) 2011 Jonathan Currie (I2C2)