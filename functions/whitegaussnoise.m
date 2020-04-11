%
% WHITEGAUSSNOISE Generate gaussian-distributed noise with uniform power
%                 spectrum distribution
%
%   x = WHITEGAUSSNOISE(N,level)
%   Generates a N-point vector (x) of Gaussian distributed random noise. The
%   standard deviation of the noise is determined by (level).
%
%   x = WHITEGAUSSNOISE(t,level)
%   The time axis (t) can be passed instead of the number of points.
%
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function noise = whitegaussnoise(N,level,matchlevel)

if nargin<3 || isempty(matchlevel)
    matchlevel = false;
end

if numel(N)>1
   N = numel(N); 
end

% Validate input
validateattributes(N,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'N')
validateattributes(level,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'level')
validateattributes(matchlevel,{'logical'},{'nonempty'},mfilename,'matchlevel')

if matchlevel
    %Generate Gaussian noise with matched standard deviation
    randvec = randn(N,1);
    %Increase amplitude of the noise vector until standard deviation matches
    %the requested noise level exactly
    noise = 0*randvec;
    amp = 0;
    while std(noise)<level
        amp = amp + level/100;
        noise = amp*randvec;
    end
else
    % Generate Gaussian noise with proper distribution standard deviation
    noise = level*randn(N,1);
end
end