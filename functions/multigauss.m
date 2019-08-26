%
% MULTIGAUSS Multigauss fitting of a distance distribution
%
%   P = MULTIGAUSS(S,K,r,Ngauss)
%   Fits the dipolar signal (S) to a distance distribution (P)using a
%   multi-gauss parametric model according to the dipolar kernel (K) and
%   distance axis (r). The function chooses the optimal number of Gaussian
%   distributions up to a maximum number given by (Ngauss) by means of the 
%   corrected Aikaike information criterion (AICC).
%
%   P = MULTIGAUSS(...,'Property',Value)
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


function [FitDistribution,FitParam,optimum,metrics] = multigauss(Signal,Kernel,DistanceAxis,maxGaussians,varargin)

%Validate user input (first three inputs are validated in lower-level functions)
if nargin<4
    error('Not enough input arguments')
end
validateattributes(maxGaussians,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'maxGaussians')

%Preallocate empty cell array for multigaussian models
multiGaussModel = cell(maxGaussians,1);

%Start with one gaussian parametric model
multiGaussModel{1} = @onegaussian;
for i=2:maxGaussians
    %And iteratively combine with additional single gaussian parametric models
   multiGaussModel{i} =  mixmodels({multiGaussModel{i-1},@onegaussian});
end

%Run optimization to see which multigauss model is the best 
[optimum,metrics] = selectmodel(multiGaussModel,Signal,DistanceAxis,Kernel,'aicc',varargin);

%Fit the data to the optimal multigauss parametric model
[FitDistribution,FitParam] = fitparamodel(Signal,Kernel,DistanceAxis,multiGaussModel{optimum},[],varargin);

return

