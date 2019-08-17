%
% WHITENOISE Generate gaussian-distributed noise with uniform power
%            spectrum distribution
%
%   x = WHITENOISE(N,level)
%   Generates a N-point vector (x) of Gaussian distributed random noise. The
%   standard deviation of the noise is determined by the (level) input
%   argument.
%
%   x = WHITENOISE(N,level,seed)
%   A different random number generator seed (default = 2) can
%   be passed as the (seed) argument.
%
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function ampnoise = whitenoise(N,level,seed)

%Validate input
if nargin<3
    seed  = 2;
end
validateattributes(N,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'seed')
validateattributes(level,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'seed')
validateattributes(seed,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'seed')

%Fix the random number generator
rng(seed);
%Generate Gaussian noise
noise = randn(N,1);
%Increase amplitude of the noise vector until standard deviation matches
%the requested noise level
ampnoise = 0;
amp = 0;
while std(ampnoise)<level
    amp = amp+0.0001;
    ampnoise = amp*noise;
end

end