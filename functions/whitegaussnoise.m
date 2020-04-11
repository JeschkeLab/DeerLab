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


function noise = whitegaussnoise(N,level,rescale)

if nargin<2 || nargin>3
    error('Two inputs (N, level) are required.');
end

if nargin<3
    rescale = '';
end

if numel(N)>1
   N = numel(N); 
end

% Validate inputs
validateattributes(N,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'N')
validateattributes(level,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'level')
validateattributes(rescale,{'char'},{},mfilename,'rescale')

rescale = strcmp(rescale,'rescale');

noise = randn(N,1);
if rescale
  noise = noise/std(noise); 
end
noise = level*noise;

end