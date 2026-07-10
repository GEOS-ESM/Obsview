function str = latstr(lat,prec)

% LATSTR Create a latitude string.

if nargin==1, prec = 4; end
if lat>=0, str = [num2str(lat,prec) 'N']; else, str = [num2str(-lat,prec) 'S']; end
