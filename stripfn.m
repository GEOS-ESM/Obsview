function [path,name,ext] = stripfn(file);

% STRIPFN Strip file name.
%
% [path, name, ext] = stripfn(filename) returns the path, name, and
%       extension for filename. For example:
%
% filename:                  path:             name:           ext:
%
% 'data/test/klklk.dat'      'data/test/'      'klklk'         '.dat'
% 'c:\data\test\klklk.dat'   'c:\data\test\'   'klklk'         '.dat'
% 'klklk.dat'                []                'klklk'         '.dat'
% 'klklk'                    []                'klklk'         []

% 10Nov97 Dick Dee

path = []; name = file; ext = [];

k = findstr(file,filesep);
if any(k),
   path = name(1:k(length(k)));
   name = name(k(length(k))+1:length(name));
end

k = findstr(name,'.');
if any(k),
   ext  = name(k(length(k)):length(name));
   name = name(1:k(length(k))-1);
end
