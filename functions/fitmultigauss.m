%
% FITMULTIGAUSS Multigauss fitting of a distance distribution
%
%   P = FITMULTIGAUSS(S,K,r,Ngauss)
%   Fits the dipolar signal (S) to a distance distribution (P)using a
%   multi-gauss parametric model according to the dipolar kernel (K) and
%   distance axis (r). The function chooses the optimal number of Gaussian
%   distributions up to a maximum number given by (Ngauss) by means of the 
%   corrected Aikaike information criterion (AICC).
%
%   [P,param,opt,metrics] = FITMULTIGAUSS(...)
%   If requested alongside the distribution (P), the optimal fit model 
%   parameters (param), the optimal number of gaussians (opt) and
%   evaluated selection metrics (metrics) are returned.
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


function [FitP,FitParam,optimum,metrics] = fitmultigauss(S,K,r,maxGaussians,method,varargin)

if nargin<5 
    method = 'aicc';
end

%Validate user input (first three inputs are validated in lower-level functions)
if nargin<4
    error('Not enough input arguments')
end
validateattributes(maxGaussians,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'maxGaussians')

%Preallocate empty cell array for multigaussian models
multiGaussModel = cell(maxGaussians,1);

%Start with one gaussian parametric model
multiGaussModel{1} = @rd_onegaussian;
for i=2:maxGaussians
    %And iteratively combine with additional single gaussian parametric models
   multiGaussModel{i} =  mixmodels({multiGaussModel{i-1},@rd_onegaussian});
end

%Run optimization to see which multigauss model is the best 
[optimum,metrics] = selectmodel(multiGaussModel,S,r,K,method,varargin);

%Fit the data to the optimal multigauss parametric model
[FitParam,FitP] = fitparamodel(S,multiGaussModel{optimum},r,K,varargin);

return

