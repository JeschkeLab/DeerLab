%
% MIXMODELS Combine parametric models into one
%
%   newmodel = MIXMODELS({@model1,@model2,...,@modelN})
%   Combines the parametric model function handles (@model1,...,@modelN)
%   into a new parametric model function handle (newmodel). The models
%   must be passed as a cell array of function handles.
%
%   The parametric model function must be of the type as the models distributed
%   in DeerLab2. The returned function handle can be used for
%   parametric model fitting as the other models.
%
%   Example: dd_mix = MIXMODELS({@dd_gauss,@dd_gauss})
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.


function finalModel = mixmodels(models)

if ~isa(models,'cell') || any(~cellfun(@(models)isa(models,'function_handle'),models))
   error('Input argument must be a cell array of valid function handles.')
end

%Check how many models are to be mixed
nModels = length(models);
% Combine the information structures of the models
%------------------------------------------------------
%Account for weight parameters between the models
mixedInfo.nparam = nModels-1;
mixedInfo.model = 'Mixed model';
for j=1:mixedInfo.nparam
    mixedInfo.parameters(j).name = sprintf('Relative amplitude of model #%i',j);
    mixedInfo.parameters(j).range = [0 1];
    mixedInfo.parameters(j).default = 1/(mixedInfo.nparam + 1);
    mixedInfo.parameters(j).units = '';
end
%Compile necessary parameters for each mode into an array
paramsplit = zeros(length(models),1);
paramsplit(1) = mixedInfo.nparam;
%Get information strucutre from each model
for i=1:length(models)
    currentModel = models{i};
    info = currentModel();
    paramsplit(i+1) =  paramsplit(i) + info.nparam;
    mixedInfo.models{i} =  info.model;
    for j=1:info.nparam
        mixedInfo.parameters(length(mixedInfo.parameters)+1).name = sprintf('%s of model #%i',info.parameters(j).name,i);
        mixedInfo.parameters(length(mixedInfo.parameters)).range = info.parameters(j).range;
        mixedInfo.parameters(length(mixedInfo.parameters)).default = info.parameters(j).default;
        mixedInfo.parameters(length(mixedInfo.parameters)).units = info.parameters(j).units;
    end
    mixedInfo.nparam =  mixedInfo.nparam + info.nparam;
end

% Merge the parametric model fuctions
%------------------------------------------------------

%Start with a dummy function handle
mixedModel = @(r,param) 0*r;

for j=1:nModels-1
    currentModel = models{j};
    mixedModel = @(r,param) (mixedModel(r,param) + param(j)*currentModel(r,param(paramsplit(j)+1:paramsplit(1+j))) );
end
%For the last model use constrained amplitude
currentModel = models{end};
mixedModel = @(r,param) (mixedModel(r,param) + max(1 - sum(param(1:paramsplit(1))),0)*currentModel(r,param(paramsplit(end-1)+1:paramsplit(end))) );

%Bundle everyhting into a function handle which can handle varargin cases
finalModel = @mixedFunction;

    %Function to allow request of information structure or model values
    function output = mixedFunction(varargin)

        if nargin==0
            output = mixedInfo;
            return
        else
            if nargin<2
                error('At least two input argumetns (ax,param) are required.')
            elseif  nargin>3
                error('Too many input argumetns given.')
            end
        
            r = varargin{1};
            if ~iscolumn(r)
               r = r.'; 
            end
            distr = mixedModel(r,varargin{2});
            output = distr;
        end
        
    end

end
