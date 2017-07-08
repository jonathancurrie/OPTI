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