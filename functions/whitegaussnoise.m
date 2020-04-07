%
% WHITEGAUSSNOISE Generate gaussian-distributed noise with uniform power
%                 spectrum distribution
%
%   x = WHITEGAUSSNOISE(N,level)
%   Generates a N-point vector (x) of Gaussian distributed random noise. The
%   standard deviation of the noise is determined by the (level) input
%   argument.
%
%   x = WHITEGAUSSNOISE(t,level)
%   The time axis (t) can be passed as well as an argument instead of the
%   number of points.
%
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function ampnoise = whitegaussnoise(N,level)

if numel(N)>1
   N = numel(N); 
end

%Validate input
validateattributes(N,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'N')
validateattributes(level,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'seed')

%Generate Gaussian noise
noise = randn(N,1);
%Increase amplitude of the noise vector until standard deviation matches
%the requested noise level
ampnoise = 0*noise;
amp = 0;
while std(ampnoise)<level
    amp = amp + level/100;
    ampnoise = amp*noise;
end

end