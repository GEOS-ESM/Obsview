function ok = isodsfile(file)

% ISODSFILE True for ODS file.
%
%    ISODSFILE(FILE) returns 1 if FILE is a readable ODS file and 0 otherwise.

% 25Mar2002 Dick Dee

ok = (exist(file,'file')==2);

if ~ok, return, end

try
    evalc(['getodsinfo ' file]);
catch
    ok = false;
end
