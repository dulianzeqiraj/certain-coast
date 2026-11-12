function tf = issueAlarm(forecast, threshold, margin)
%ISSUEALARM Alarm rule: true where the forecast plus the margin reaches the threshold.
%   TF = certaincoast.alarms.issueAlarm(FORECAST, THRESHOLD, MARGIN) returns a logical array of the size of FORECAST.
%
%   Example (chloride, mg/L)
%       alarm = certaincoast.alarms.issueAlarm(newForecast, 250, m);
%
%   See also certaincoast.alarms.conformalMargin.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

tf = double(forecast) + margin >= threshold;
end
