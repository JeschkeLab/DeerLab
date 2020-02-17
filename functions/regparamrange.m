%
% REGPARAMRANGE Regularization parameter range estimator
%
%   alphas = REGPARAMRANGE(K,L)
%   Estimates an array of regularization parameter candidates (alphas) from
%   the generalized singular value decomposition (GSVD) of the experiment
%   kernel matrix (K) and regularization operator matrix (L).
%
%   alphas = REGPARAMRANGE(...,'Property',Value)
%   Additional (optional) arguments can be passed as property-value pairs.
%
% The property-value pairs to be passed as options can be set in any order. 
%
%   'Resolution' - Resolution of the array of alpha candidates, on a 
%                  base-10 logarithmic scale (default 0.1, corresponding to
%                  a tenth of a decade).
%
%   'NoiseLevel' - Estimate of the noise standard deviation of the
%                  signal; used for scaling of the singular values (default=0).
%
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function alphas = regparamrange(K,L,varargin)

% Check if user requested some options via name-value input
[stdNoise,lgResolution] = parseoptional({'NoiseLevel','Resolution'},varargin);

if isempty(stdNoise)
    stdNoise = 0;
end
if isempty(lgResolution)
    lgResolution = 0.1;  % resolution of log10(alpha)
end
if issparse(L)
    L = full(L); 
end
validateattributes(K,{'numeric'},{'nonempty'})
validateattributes(L,{'numeric'},{'nonempty'})
validateattributes(lgResolution,{'numeric'},{'scalar'})
validateattributes(stdNoise,{'numeric'},{'scalar'})

% Set alpha range
%-------------------------------------------------------------------------------
minmax_ratio = 16*eps*1e6;  % ratio of smallest to largest alpha

% Scaling by noise. This improves L curve corner detection for DEER.
minmax_ratio = minmax_ratio*2^(stdNoise/0.0025);

% Get generalized singular values of K and L
%-------------------------------------------------------------------------------
singularValues = gsvd(K,L,0);
DerivativeOrder = size(L,2) - size(L,1); % get order of derivative (=number of inf in singval)
singularValues = singularValues(1:end-DerivativeOrder); % remove inf 
singularValues = singularValues(end:-1:1); % sort in decreasing order
lgsingularValues = log10(singularValues);

% Calculate range based on singular values
%-------------------------------------------------------------------------------
lgrangeMax = lgsingularValues(1);
lgrangeMin = max([lgsingularValues(end),lgsingularValues(1)+log10(minmax_ratio)]);
lgrangeMax = floor(lgrangeMax/lgResolution)*lgResolution;
lgrangeMin = ceil(lgrangeMin/lgResolution)*lgResolution;
if lgrangeMax < lgrangeMin
    warning('Determined range maximum is smaller than range minimum.');
    temp = lgrangeMax;
    lgrangeMax = lgrangeMin;
    lgrangeMin = temp;
end
lgalpha = lgrangeMax:-lgResolution:lgrangeMin;
alphas = 10.^lgalpha;
alphas = alphas(:);
