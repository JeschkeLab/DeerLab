%
% CORRECTZEROTIME Zero-time correction of dipolar spectroscopy signals
%
%   tc = CORRECTZEROTIME(V,t)
%   Determines the zero time of a dipolar signal (V) and corrects
%   the time axis (t) for it, returning a corrected time axis (tc).
%   If t is in ns/us, tc will be be in ns/us as well.
%
%   tc = CORRECTZEROTIME(V,t,t0)
%   Corrects the time axis (t) for a given zero-time (t0), returning the
%   corrected time axis.
%
%   [tc,t0,pos] = CORRECTZEROTIME(V,t)
%   Returns the corrected time axis (tc), the zero-time (t0) in ns/us and the 
%   zero-time array index (pos). 
%
%

% This file is a part of DeerAnalysis. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.


function [tcorr,t0,idxt0] = correctzerotime(V,t,t0)

if nargin<2
    error('Not enough input arguments.')
end

if nargin<3 || isempty(t0)
    t0 = [];
else
    validateattributes(t0,{'numeric'},{'scalar','nonnegative'},mfilename,'ZeroTime')
end
if ~iscolumn(V)
    V = V.';
end
if any(t<0)
   t = t + abs(min(t)); 
end
validateattributes(V,{'numeric'},{'2d'},mfilename,'S')
validateattributes(t,{'numeric'},{'nonempty'},mfilename,'t')

%Generate finely-grained interpolated signal and time axis
resolution = 4;
tfine = linspace(min(t),max(t),(numel(t)-1)*resolution+1);
Vfine = interp1(t,real(V),tfine,'spline');
% get zero time, if not provided
if isempty(t0)
    % Determine maximum
    [~,maxPos] = max(Vfine);
    idxt0 = 1;
    %If maximum is not the first or last point, then do moment analysis
    if maxPos>1 && maxPos<length(Vfine)
        % Determine the width of the interval to use for the integral
        if maxPos<length(Vfine)/2
          maxDelta = floor((maxPos-1)/2);
        else
          maxDelta = floor((length(Vfine)-maxPos)/2);
        end
        offsetrange = -maxDelta:maxDelta;
        
        %Look around the maximum for lowest moment integral
        minIntegral = realmax;
        idxt0 = maxPos;
        for idx = maxPos-maxDelta:maxPos+maxDelta
            %Query new candidate for zero position and integrate
            dt = tfine(idx+offsetrange)-tfine(idx);
            Integral = sum(Vfine(idx+offsetrange).*dt);
            %If integral is lower than prior best, then update new candidate
            if abs(Integral) < minIntegral
                minIntegral = abs(Integral);
                idxt0 = idx;
            end
        end
    end
    t0 = tfine(idxt0);
end

% Correct time axis
tcorr = t - t0;
[~,idxt0] = min(abs(tcorr));

end
