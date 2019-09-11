%
% FITMULTIGAUSS Multi-Gauss fitting of a distance distribution
%
%   P = FITMULTIGAUSS(S,K,r,Ngauss,method)
%   Fits a multi-Gauss parametric distance distribution model to the dipolar
%   signal (S), using the dipolar kernel (K) and distance axis (r). The function
%   compares multi-Gaussian distributions with up to a maximum number of Gaussians
%   given by (Ngauss) and determines the optimum one using the model selection
%   criterion given in (method) ('AIC', 'BIC', or 'AICc'). The fitted
%   distribution is returned in P.
%
%   [P,param,opt,metrics,Peval] = FITMULTIGAUSS(...)
%   If requested alongside the distribution (P), the optimal fit model 
%   parameters (param), the optimal number of Gaussians (opt) and
%   evaluated selection metrics (metrics) are returned. The fitted distance
%   distributions fitted for ech multigauss model can be requested as a
%   fifth output argument (Peval).
%
%   P = FITMULTIGAUSS(...,'Property',Value)
%   Additional (optional) arguments can be passed as property-value pairs. 
%   
%   See "help fitparamodel" for a detailed list of the property-value pairs
%   accepted by the function.
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function [Pfit,param,nGaussOpt,metrics,Peval] = fitmultigauss(S,K,r,maxGaussians,method,varargin)

% Validate user input (S, K, r, and method are validated in lower-level functions)
if nargin<4
    error('Not enough input arguments.')
else
    validateattributes(maxGaussians,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'maxGaussians')
end
if nargin<5 
    method = 'aicc';
end

% Compile list of multi-Gaussian models
multiGaussModels = cell(maxGaussians,1);
multiGaussModels{1} = @rd_onegaussian;
if maxGaussians>=2, multiGaussModels{2} = @rd_twogaussian; end
if maxGaussians>=3, multiGaussModels{3} = @rd_threegaussian; end
if maxGaussians>=4, multiGaussModels{4} = @rd_fourgaussian; end
if maxGaussians>=5, multiGaussModels{5} = @rd_fivegaussian; end
for i = 6:maxGaussians
   multiGaussModels{i} =  mixmodels({multiGaussModels{i-1},@rd_onegaussian});
end

% Run fitting and model selection to see which multi-Gauss model is optimal
[nGaussOpt,metrics,fitparams] = selectmodel(multiGaussModels,S,r,K,method,varargin);

% Calculate the distance distribution for the optimal multi-Gauss model
param = fitparams{nGaussOpt};
optModel = multiGaussModels{nGaussOpt};
Pfit = optModel(r,param);

if nargout>4
    for i=1:maxGaussians
        Peval(i,:) = multiGaussModels{i}(r,fitparams{i});
    end
end

return
