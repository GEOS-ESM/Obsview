function varargout = dconfig(varargin)

% DCONFIG Load system-specific parameters.

gsiinfo   % load descriptions of data attributes etc.

nargin = length(varargin);

if nargin==0, whos; end

for i = 1:length(varargin),

  varname = char(varargin{i});
  eval(['varargout{i} = ' varname ';'])
  
end

