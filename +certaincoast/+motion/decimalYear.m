function t = decimalYear(y, m, d)
%DECIMALYEAR Date as a decimal year: the middle of the day, as a fraction of its own year.
%   T = certaincoast.motion.decimalYear(Y, M, D) for arrays of years, months and days of the same size (or scalars)
%   returns Y + (day of the year - 1/2) / (days in the year), so 2019-11-26 is 2019.9027 and leap years count 366 days.
%   T = certaincoast.motion.decimalYear(YYYYMMDD) accepts the dates as numbers such as 20191126.
%
%   Example
%       t = certaincoast.motion.decimalYear(2019, 11, 26);     % the Durres earthquake of 26 November 2019
%
%   See also certaincoast.motion.fitMotionSeries, certaincoast.motion.readSeriesCsv.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin == 1
    v = y;
    y = floor(v / 10000);
    m = floor(mod(v, 10000) / 100);
    d = mod(v, 100);
end
if isscalar(y), y = y + zeros(size(m)); end
if isscalar(m), m = m + zeros(size(y)); end
if isscalar(d), d = d + zeros(size(y)); end
if ~isequal(size(y), size(m), size(d)) || any(m(:) < 1 | m(:) > 12 | d(:) < 1 | d(:) > 31)
    error('certaincoast:motion:badDate', 'need years, months 1-12 and days 1-31 of the same size.');
end
first = datenum(y, 1, 1);
len = datenum(y + 1, 1, 1) - first;
t = y + (datenum(y, m, d) - first + 0.5) ./ len;
end
