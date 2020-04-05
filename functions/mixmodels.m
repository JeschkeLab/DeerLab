%
% MIXMODELS Combine parametric models into one
%
%   newmodel = MIXMODELS(@model1,@model2,...,@modelN)
%   Combines the parametric model function handles (@model1,...,@modelN)
%   into a new parametric model function handle (newmodel). The models
%   must be passed as a cell array of function handles.
%
%   The parametric model function must be of the type as the models distributed
%   in DeerLab2. The returned function handle can be used for
%   parametric model fitting as the other models.
%
%   Example: dd_mix = MIXMODELS(@dd_gauss,@dd_gauss)
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function mixModelFcn = mixmodels(varargin)

if numel(varargin)==1
    models = {varargin};
else
    models = varargin;
end

if numel(varargin)==0
    error('At least one model must be provided.')
end

if ~all(cellfun(@(M)isa(M,'function_handle'),models))
    error('Input arguments must all be function handles.')
end

% Detemine number of models to be mixed
nModels = numel(models);

% Combine the information structures of the models
%-------------------------------------------------------------------------------
Info.model = 'Mixed model';

% Add amplitudes for each model except last
for j = 1:nModels-1
    Info.parameters(j).name = sprintf('Relative amplitude of model #%i',j);
    Info.parameters(j).range = [0 1];
    Info.parameters(j).default = 1/nModels;
    Info.parameters(j).units = '';
end
Info.nparam = nModels-1; % amplitudes

% Combine parameters for each model into an array
paramsplit = zeros(numel(models),1);
paramsplit(1) = Info.nparam;

% Add information structure from each model
for i = 1:length(models)
    info = models{i}();
    paramsplit(i+1) =  paramsplit(i) + info.nparam;
    Info.models{i} = info.model;
    for j = 1:info.nparam
        Info.parameters(length(Info.parameters)+1).name = sprintf('%s of model #%i',info.parameters(j).name,i);
        Info.parameters(length(Info.parameters)).range = info.parameters(j).range;
        Info.parameters(length(Info.parameters)).default = info.parameters(j).default;
        Info.parameters(length(Info.parameters)).units = info.parameters(j).units;
    end
    Info.nparam =  Info.nparam + info.nparam;
end

% Merge the parametric model functions
%-------------------------------------------------------------------------------

% Start with a dummy function handle
mixedModel = @(r,param) 0*r;
% Add models one by one using nested function handles
for j = 1:nModels-1
    currentModel = models{j};
    mixedModel = @(r,param) (mixedModel(r,param) + param(j)*currentModel(r,param(paramsplit(j)+1:paramsplit(1+j))) );
end
% For the last model use constrained amplitude
currentModel = models{end};
mixedModel = @(r,param) (mixedModel(r,param) + max(1 - sum(param(1:paramsplit(1))),0)*currentModel(r,param(paramsplit(end-1)+1:paramsplit(end))) );

% Bundle everyhting into a function handle which can handle varargin cases
mixModelFcn = @mixedFunction;

    % Function to allow request of information structure or model values
    function output = mixedFunction(varargin)

        if nargin==0
            output = Info;
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
