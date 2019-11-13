%
% SPECUNITS Specify units for time/distance axis vectors
%
%   t = SPECUNITS(t,'ns')
%   t = SPECUNITS(t,'us')
%   Checks the time-axis (t) and converts it to ensure other functions
%   recognize it with the specified unit ('ns' or 'us').
%
%   t = SPECUNITS(t,'nm')
%   t = SPECUNITS(t,'A')
%   Checks the distance-axis (r) and converts it to ensure other functions
%   recognize it with the specified unit ('nm' or 'A').
%

% This file is a part of DeerAnalysis. License is MIT (see LICENSE.md).
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function [ax,unit] = specunits(ax,type)

%Check the type of unit being requested by user
switch type
    case 'ns'
        tmax = max(ax) - min(ax);
        if tmax<=30
            ax = ax/1000;
            unit = 'us';
        else
            unit = 'ns';
        end
        
    case 'us'
        tmax = max(ax) - min(ax);
        if tmax>30
            ax = ax*1000;
            unit = 'ns';
        else
            unit = 'us';
        end
        
    case 'nm'
        if max(ax)>=20
            ax = ax*10;
            unit = 'A';
        else
            unit = 'nm';
        end
        
    case 'A'
        if max(ax)<20
            ax = ax/10;
            unit = 'nm';
        else
            unit = 'A';
        end
        
end

%Correct for round-off errors
ax = round(ax,14);


end