function [correctedSignal,correctedTimeAxis,ZeroTime] = correctZeroTime(Signal,TimeAxis,ZeroTime)

%Generate finely-grained interpolated signal and time axis
FineTimeAxis = min(TimeAxis):1:max(TimeAxis);
FineSignal=interp1(TimeAxis,real(Signal),FineTimeAxis,'spline',real(Signal(1)));
% get zero time, if not provided
FineSignal=real(FineSignal);
if nargin<3
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
          SmallestDistance=MaxEndDistance;
        else
          SmallestDistance=NewPosDistance;
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
    ZeroTime=FineTimeAxis(ZeroTimePos);
else
    [~,ZeroTimePos] = min(abs(FineTimeAxis-ZeroTime));
end

% Get a smoothed estimate of the zero-position
if ZeroTimePos > 10
  SmoothingAxis = 1:min(2*ZeroTimePos,length(FineSignal));
  [p,~] = polyfit(SmoothingAxis,FineSignal(SmoothingAxis),11);
  SmoothZeroTimeRegion = polyval(p,SmoothingAxis);
  correctedSignal = Signal/SmoothZeroTimeRegion(ZeroTimePos);
else
  correctedSignal = Signal;
end

% Correct time axis
correctedTimeAxis = TimeAxis - ZeroTime*ones(size(TimeAxis));

