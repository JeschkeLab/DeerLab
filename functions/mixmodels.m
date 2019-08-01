function finalModel = mixmodels(models)

if ~isa(models,'cell') || any(~cellfun(@(models)isa(models,'function_handle'),models))
   error('Input argument must be a cell array of valid function handles.')
end

%Check how many models are to be mixed
nModels = length(models);
% Combine the information structures of the models
%------------------------------------------------------
%Account for weight parameters between the models
mixedInfo.nParam = nModels-1;
mixedInfo.Model = 'Mixed model';
for j=1:mixedInfo.nParam
    mixedInfo.parameters(j).name = sprintf('Relative amplitude of model #%i',j);
    mixedInfo.parameters(j).range = [0 1];
    mixedInfo.parameters(j).default = 1/(mixedInfo.nParam + 1);
    mixedInfo.parameters(j).units = '';
end
%Compile necessary parameters for each mode into an array
paramsplit = zeros(length(models),1);
paramsplit(1) = mixedInfo.nParam;
%Get information strucutre from each model
for i=1:length(models)
    currentModel = models{i};
    info = currentModel();
    paramsplit(i+1) =  paramsplit(i) + info.nParam;
    mixedInfo.models{i} =  info.Model;
    mixedInfo.equations{i} =  info.Equation;
    for j=1:info.nParam
        mixedInfo.parameters(length(mixedInfo.parameters)+1).name = sprintf('%s of model #%i',info.parameters(j).name,i);
        mixedInfo.parameters(length(mixedInfo.parameters)).range = info.parameters(j).range;
        mixedInfo.parameters(length(mixedInfo.parameters)).default = info.parameters(j).default;
        mixedInfo.parameters(length(mixedInfo.parameters)).units = info.parameters(j).units;
    end
    mixedInfo.nParam =  mixedInfo.nParam + info.nParam;
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
            distr = mixedModel(varargin{1},varargin{2});
            output = distr;
        end
        
    end

end
