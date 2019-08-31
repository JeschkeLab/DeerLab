%
% CORRECTZEROTIME Zero-time correction of dipolar spectroscopy signals
%
%   tc = CORRECTZEROTIME(V,t)
%   Determines the zero time of a dipolar signal (V) and corrects
%   the time axis (t) for it, returning a corrected time axis (tc).
%   If t is in ns/us, tc will be be in ns/us as well.
%
%   tc = CORRECTZEROTIME(V,t,zt)
%   Corrects the time axis (t) for a given zero-time (zt), returning the
%   corrected time axis.
%
%   [tc,zt,pos] = CORRECTZEROTIME(V,t)
%   Returns the corrected time axis (tc), the zero-time (zt) in ns/us and the 
%   zero-time array index (pos). 
%
% Adapted from Gunnar Jeschke, DeerAnalysis 2018
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function [correctedTimeAxis,ZeroTime,ZeroTimePos] = correctzerotime(Signal,TimeAxis,ZeroTime)

if nargin<2
    error('Not enough input arguments.')
end

if nargin<3 || isempty(ZeroTime)
    ZeroTime = [];
else
    validateattributes(ZeroTime,{'numeric'},{'scalar','nonnegative'},mfilename,'ZeroTime')
end
if ~iscolumn(Signal)
    Signal = Signal.';
end
validateattributes(Signal,{'numeric'},{'2d'},mfilename,'Signal')
validateattributes(TimeAxis,{'numeric'},{'nonnegative','nonempty'},mfilename,'TimeAxis')

usesMicrosecondUnits = mean(abs(diff(TimeAxis))) < 1;

if usesMicrosecondUnits
    TimeAxis = TimeAxis*1000; % convert us -> ns
end

%Generate finely-grained interpolated signal and time axis
FineTimeAxis = min(TimeAxis):1:max(TimeAxis);
FineSignal=interp1(TimeAxis,real(Signal),FineTimeAxis,'spline',real(Signal(1)));
% get zero time, if not provided
FineSignal=real(FineSignal);
if isempty(ZeroTime)
    % Determine maximum
    [~,maxPos]=max(FineSignal);
    ZeroTimePos=1;
    %If maximum is not the first point in signal, then do moment-analysis
    if maxPos>1 && maxPos<length(FineSignal)
        %
        % Procedure schematic:
        %
        %   NewPosDistance                 MaxEndDistance
        %  <------------->     <--------------------------------------->
        % 1            newPos maxPos                                  end
        % |--------------|-----|---------------------------------------|
        %         <--|--------------------> 2xMaxAllowedDistance
        %         TrialPos
        % <-----------------------> 2xMaxAllowedDistance
        %        Integration
        %
        %
        %Determine the maximum allowed distance for search
        NewPos = maxPos - 1;
        MaxEndDistance = length(FineSignal) - maxPos;
        NewPosDistance = NewPos;
        if MaxEndDistance<NewPosDistance
            SmallestDistance = MaxEndDistance;
        else
            SmallestDistance = NewPosDistance;
        end
        MaxAllowedDistance = floor(NewPosDistance/2);
        BestIntegral = 1e20;
        %Look around the maximum for lowest signal integral over MaxAllowedDistance
        ZeroTimePos = maxPos;
        for TrialPos = -MaxAllowedDistance+maxPos:MaxAllowedDistance+maxPos
            %Query new candidate for zero position and integrate
            Integral = 0;
            for j = -MaxAllowedDistance:MaxAllowedDistance
                Integral = Integral + FineSignal(TrialPos+j)*j;
            end
            %If integral is lower than prior best, then update new candidate
            if abs(Integral) < BestIntegral
                BestIntegral = abs(Integral);
                ZeroTimePos = TrialPos;
            end
        end
    end
    ZeroTime = FineTimeAxis(ZeroTimePos);
else
    [~,ZeroTimePos] = min(abs(FineTimeAxis-ZeroTime));
end

% Correct time axis
correctedTimeAxis = TimeAxis - ZeroTime;

if usesMicrosecondUnits
    ZeroTime = ZeroTime/1000; % convert ns -> us
    correctedTimeAxis = correctedTimeAxis/1000;  % convert ns -> us
end

end
