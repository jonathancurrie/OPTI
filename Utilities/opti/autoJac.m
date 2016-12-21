function jac = autoJac(fun,x0,varargin)
%AUTOJAC  Returns a Automatically Differentiated (AD) Jacobian
% 
%   jac = autoJac(fun,x0) uses the user supplied Matlab function to
%   generate the Jacobian using automatic differentiation. The user
%   function must be a Matlab function and must not call external code or
%   class / toolbox functions. If your function breaks any of these
%   conditions consider using mklJac instead.
%
%   The underlying AD algorithm is adiff by William McIlhagga and its
%   documentation pdf is provided in the Differentiation folder. See the
%   BSD license below the code.

%   Copyright (C) 2011 Jonathan Currie (IPL)

if(~isa(fun,'function_handle'))
    error('Fun should be a function handle!');
end

[~,jac] = adiffget(fun(adiff(x0),varargin{:}));

end

% Copyright (c) 2010, William McIlhagga
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
