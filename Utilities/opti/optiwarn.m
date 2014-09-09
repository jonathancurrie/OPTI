function optiwarn(identifier,message,varargin)
%Display a warning without the backtrace

s = warning('backtrace');
warning('off','backtrace');
warning(identifier,message,varargin{:});
warning(s);