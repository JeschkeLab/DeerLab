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

