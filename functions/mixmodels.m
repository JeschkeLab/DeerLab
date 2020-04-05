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
Info.models = {};

% Add amplitudes for each model except last
for j = 1:nModels-1
    Info.parameters(j).name = sprintf('Component %i: Relative amplitude',j);
    Info.parameters(j).range = [0 1];
    Info.parameters(j).default = 1/nModels;
    Info.parameters(j).units = '';
end
pidx_amp = 1:nModels-1;

% Combine info structures from all models
idx = pidx_amp(end);
for i = 1:nModels
    info = models{i}();
    pidx{i} = idx + (1:info.nparam);
    idx = idx + info.nparam;
    Info.models{i} = info.model;
    param_ = info.parameters;
    for j = 1:info.nparam
        param_(j).name = sprintf('Component %i: %s',i,param_(j).name);
    end
    Info.parameters = [Info.parameters param_];
end

Info.nparam = numel(Info.parameters);

% Mixed model function handle
%-------------------------------------------------------------------------------
mixModelFcn = @mixedFunction;

    % Function to allow request of information structure or model values
  function output = mixedFunction(varargin)

        if nargin==0
            output = Info;
            return
        end
        
        if nargin<2
            error('At least two input arguments (r,param) are required.')
        elseif  nargin>3
            error('Only two input arguments are allows.')
        end
        
        x = varargin{1};
        params = varargin{2};
        if ~iscolumn(x)
            x = x.';
        end
        
        amp = params(pidx_amp);
        amp(end+1) = max(1-sum(amp),0);
        
        y = 0;
        for k = 1:numel(models)
            y = y + amp(k)*models{k}(x,params(pidx{k}));
        end
        output = y;
        
    end
    
end
