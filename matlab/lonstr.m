function str = lonstr(lon,prec)

% LONSTR Create a longitude string.

if nargin==1, prec = 4; end
if lon<=0, str = [num2str(-lon,prec) 'W']; else, str = [num2str(lon,prec) 'E']; end
