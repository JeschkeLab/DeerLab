%
% WHITEGAUSSNOISE Generate gaussian-distributed noise with uniform power
%                 spectrum distribution
%
%   x = WHITEGAUSSNOISE(N,level)
%   Generates a N-point vector (x) of Gaussian distributed random noise. The
%   standard deviation of the noise is determined by the (level) input
%   argument.
%
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function ampnoise = whitegaussnoise(N,level)

%Validate input
validateattributes(N,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'seed')
validateattributes(level,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'seed')

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