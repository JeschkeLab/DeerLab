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


function noise = whitegaussnoise(N,level)

if numel(N)>1
   N = numel(N); 
end

%Validate input
validateattributes(N,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'N')
validateattributes(level,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'seed')

%Generate Gaussian noise
randvec = randn(N,1);
%Increase amplitude of the noise vector until standard deviation matches
%the requested noise level
noise = 0*randvec;
amp = 0;
while std(noise)<level
    amp = amp + level/100;
    noise = amp*randvec;
end

end