%
%   ISNANOSECOND
%
%   logical = ISNANOMETER(r)
%   If the time axis is in nanseconds, the function returns (true), if
%   it is in microseconds it returns (false).
%

% This file is a part of DeerAnalysis. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function logical = isnanosecond(t)

if iscolumn(t)
    t = t.';
end

if max(t) - min(t) > 30
    logical = true;
else
    logical = false;
end

end