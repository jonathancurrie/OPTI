% COINW  Write a Mathematical File from Matlab Values
%
% THIS IS A LOW LEVEL FUNCTION - USE coinWrite() INSTEAD!
%
% coinW uses the COIN Utilities library to open and write the file.
%
%   coinW(prob,path,filetype)
%
%   Input arguments:
%       prob - problem structure (described below)
%       path - full path to the file on your PC.
%       filetype - one of the following types: [mps,lp]
%
%   Problem Structure Fields (all required):
%       H - quadratic objective matrix (sparse, may be empty)
%       f - linear objective vector
%       A - linear inequality and equality constraint matrix (sparse)
%       rl - linear constraint lower bounds
%       ru - linear constraint upper bounds
%       lb - decision variable lower bounds
%       ub - decision variable upper bounds
%       int - integer variable vector (8 bit integer / char)
%       name - problem name (string)
%       sos_type - SOS types (char array)
%       sos_index - SOS indicies (cell array of double vectors)
%       sos_weight - SOS weights (cell array of double vectors)
%       objbias - Objective bias term (optional, negated internally)
%
%
%   Copyright (C) 2011 Jonathan Currie (I2C2)

