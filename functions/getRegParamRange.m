function alpha = getRegParamRange(Kernel,RegMatrix,varargin)

%Check if user requested some options via name-value input
[NoiseDeviation,logResolution] = parseOptional({'NoiseDeviation','logResolution'},varargin);


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
