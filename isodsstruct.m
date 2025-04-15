function ok = isodsstruct(s)

% ISODSSTRUCT True for ODS structure.
%
%    ISODSSTRUCT(S) returns 1 if S is a ODS structure and 0 otherwise.
%
% See also: ODSCREATE

% 03Dec99 Dick Dee

ok = isstruct(s);

if ~ok, return, end

fields = {'first_julian_day';'latest_julian_day';'latest_synoptic_hour'};

for i = 1:length(fields),

    ok = isfield(s,fields{i});
    if ~ok, return, end

end

