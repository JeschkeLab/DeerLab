%
% REGPARAMRANGE Regularization parameter range estimator
%
%   alphas = REGPARAMRANGE(K,L)
%   Estimates an array of regularization parameter candidates (alphas) from
%   the generalized singular value decomposition (GSVD) of the dipolar
%   kernel (K) and regularization operator (L).
%
%   alphas = REGPARAMRANGE(...,'Property',Value)
%   Additional (optional) arguments can be passed as property-value pairs.
%
% The property-value pairs to be passed as options can be set in any order. 
%
%   'NoiseDeviation' - Estimation of the noise standard deviation of the
%                      signal fro scaling of the singular values (default=0).
%
%   'logResolution' - Logarithmic scale resolution of the array of alpha
%                     candidates(default=0.1).
%   
% Adapted from Stefan Stoll
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.


function alpha = regparamrange(Kernel,RegMatrix,varargin)

%Check if user requested some options via name-value input
[NoiseDeviation,logResolution] = parseoptional({'NoiseDeviation','logResolution'},varargin);


if isempty(NoiseDeviation)
  NoiseDeviation = 0;
end
if isempty(logResolution)
    logResolution = 0.1;  % resolution of log10(alpha)
end
if issparse(RegMatrix)
   RegMatrix = full(RegMatrix); 
end
validateattributes(Kernel,{'numeric'},{'nonempty'})
validateattributes(RegMatrix,{'numeric'},{'nonempty'})
validateattributes(logResolution,{'numeric'},{'scalar'})
validateattributes(NoiseDeviation,{'numeric'},{'scalar'})

% Set alpha range
%-------------------------------------------------------------
minmax_ratio = 16*eps*1e6;  % ratio of smallest to largest alpha

% Scaling by noise. This improves L curve corner detection for DEER
minmax_ratio = minmax_ratio*2^(NoiseDeviation/0.0025);

% Get generalized singular values of K and L
%-------------------------------------------------------------
singularValues = gsvd(Kernel,RegMatrix,0);
DerivativeOrder = size(RegMatrix,2) - size(RegMatrix,1); % get order of derivative (=number of inf in singval)
singularValues = singularValues(1:end - DerivativeOrder); % remove inf 

singularValues = singularValues(end:-1:1); % sort in decreasing order

% Calculate range based on singular values
%-------------------------------------------------------------
logRegParamMax = log10(singularValues(1));
logRegParamMin = log10(max([singularValues(end),singularValues(1)*minmax_ratio]));
logRegParamMax = floor(logRegParamMax/logResolution)*logResolution;
logRegParamMin = ceil(logRegParamMin/logResolution)*logResolution;
if logRegParamMax < logRegParamMin
  temp=logRegParamMax;
  logRegParamMax=logRegParamMin;
  logRegParamMin=temp;
end
lgalpha = logRegParamMax:-logResolution:logRegParamMin;
alpha = 10.^lgalpha;
