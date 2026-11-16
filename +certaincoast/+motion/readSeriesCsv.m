function S = readSeriesCsv(file)
%READSERIESCSV Read InSAR or GNSS displacement time series from a CSV file with one point per row.
%   S = certaincoast.motion.readSeriesCsv(FILE) reads a comma-separated file whose header names the columns. Every
%   column whose name is a date YYYYMMDD holds the displacement of each point on that date [mm]; this is the layout
%   of the public products of the European Ground Motion Service (EGMS) and of time series exported from LiCSBAS.
%   The other columns (point identifier, coordinates, line-of-sight unit vector, mean velocity and so on) are kept
%   as attributes. S has the fields
%       t        1-by-m decimal years of the dates (certaincoast.motion.decimalYear)
%       dates    1-by-m dates as numbers YYYYMMDD
%       D        n-by-m displacements of the n points
%       attr     struct with one field per other column: a numeric column if every value is a number, otherwise a
%                cell array of strings (field names made valid for MATLAB)
%   The reader uses core functions only; for very large products, select the area of interest before reading.
%   Public products are downloaded by the user from their portals; the repository holds no field data.
%
%   Example
%       S = certaincoast.motion.readSeriesCsv('points.csv');        % a product downloaded by the user
%       F = certaincoast.motion.fitMotionSeries(S.t, S.D(1, :));
%
%   See also certaincoast.motion.fitMotionSeries, certaincoast.motion.losToVertical.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

fid = fopen(file, 'r');
if fid < 0
    error('certaincoast:motion:file', 'cannot open %s.', file);
end
cleaner = onCleanup(@() fclose(fid));
header = fgetl(fid);
if ~ischar(header)
    error('certaincoast:motion:empty', '%s is empty.', file);
end
names = strtrim(strsplit(strtrim(header), ','));
names = regexprep(names, '^"|"$', '');
isDate = ~cellfun(@isempty, regexp(names, '^\d{8}$', 'once'));
if ~any(isDate)
    error('certaincoast:motion:noDates', 'no column of %s is named as a date YYYYMMDD.', file);
end
rows = {};
while true
    line = fgetl(fid);
    if ~ischar(line)
        break
    end
    if isempty(strtrim(line))
        continue
    end
    rows{end + 1, 1} = strsplit(line, ','); %#ok<AGROW>
end
n = numel(rows);
m = numel(names);
cells = cell(n, m);
for i = 1:n
    r = rows{i};
    if numel(r) ~= m
        error('certaincoast:motion:columns', 'row %d of %s has %d values for %d columns.', i, file, numel(r), m);
    end
    cells(i, :) = regexprep(strtrim(r), '^"|"$', '');
end
dates = cellfun(@str2double, names(isDate));
S.dates = dates;
S.t = certaincoast.motion.decimalYear(dates);
S.D = reshape(str2double(cells(:, isDate)), n, nnz(isDate));
S.attr = struct();
other = find(~isDate);
for j = other
    values = cells(:, j);
    numbers = str2double(values);
    field = matlab_name(names{j});
    if all(~isnan(numbers) | cellfun(@isempty, values) | strcmpi(values, 'nan'))
        S.attr.(field) = numbers;
    else
        S.attr.(field) = values;
    end
end
end

function f = matlab_name(s)
% a valid field name from a column name
f = regexprep(s, '[^A-Za-z0-9_]', '_');
if isempty(f) || ~isletter(f(1))
    f = ['c_' f];
end
f = f(1:min(end, 63));
end
