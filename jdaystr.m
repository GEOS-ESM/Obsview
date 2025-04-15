function str = jdaystr(jday,dateform)

% JDAYSTR Convert Julian day to calendar date.
%
%	The jdaystr function converts Julian day numbers to date strings.
%	Calling syntax is identical to DATESTR.

% 03Nov97 Dick Dee

if nargin==1,

   str = datestr(jday + datenum('May 23 1968') - 2440000);

else

   str = datestr(jday + datenum('May 23 1968') - 2440000, dateform);

end
