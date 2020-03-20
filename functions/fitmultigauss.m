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
%   P = FITMULTIGAUSS(S,t,r,Ngauss,method)
%   If a the default kernel is to be used, the time axis (t) can be passed
%   instead of the kernel.
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
%   'BckgModel' - Function handle to the background model to be fitted
%                 along the multigauss distance distribution model.
%                 Requires the time-axis to be passed instead of the
%                 kernel.
%
%   'Lower' - Array [<r>_min FWHM_min] containing the lower bound for the
%             FWHM and mean distance of all the Gaussians.
%
%   'Upper' -  Array [<r>_max FWHM_max] containing the upper bound values
%              for the FWHM and mean distance of all the Gaussians.
%
%   See "help fitparamodel" for a detailed list of other property-value pairs
%   accepted by the function.
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function [Pfit,param,nGaussOpt,metrics,Peval] = fitmultigauss(S,K,r,maxGaussians,method,varargin)


if ~license('test','optimization_toolbox')
   error('DeerLab could not find a valid licence for the Optimization Toolbox. Please install the add-on to use fitmultigauss.')
end

% Validate user input (S, K, r, and method are validated in lower-level functions)
if nargin<4
    error('Not enough input arguments.')
else
    validateattributes(maxGaussians,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'maxGaussians')
end
if nargin<5
    method = 'aicc';
end
if ~all(size(K) > 1)
    t = K;
    K = dipolarkernel(t,r);
end
warning('off','DeerLab:parseoptional')
%Parse the optional parameters in the varargin
[Upper,Lower,BckgModel] = parseoptional({'Upper','Lower','Background','internal::parseLater'},varargin);
warning('on','DeerLab:parseoptional')
if ~isempty(Upper) && isempty(BckgModel) && length(Upper)~=2
    error('Upper property must be an array [<r>_max FWHM_max]')
elseif ~isempty(Upper) && ~isempty(BckgModel) && length(Upper)<4
    error('Upper property must be an array [<r>_max FWHM_max lambda_max Bparam_max]')
end
if ~isempty(Lower) && isempty(BckgModel)  && length(Lower)~=2
    error('Lower property must be an array [<r>_min FWHM_min]')
elseif ~isempty(Upper) && ~isempty(BckgModel) && length(Lower)<4
    error('Lower property must be an array [<r>_min FWHM_min lambda_min Bparam_min]')
end
if ~isempty(BckgModel) && ~exist('t','var')
    error('Time axis must be provided for a time-domain fit.')
end
%Remove the Lower and Upper options from varargin so they are not passed to fitparamodel
Idx = find(cellfun(@(x)(ischar(x) && strcmpi(x,'upper')),varargin));
varargin(Idx:Idx+1) = [];
Idx = find(cellfun(@(x)(ischar(x) && strcmpi(x,'lower')),varargin));
varargin(Idx:Idx+1) = [];
Idx = find(cellfun(@(x)(ischar(x) && strcmpi(x,'background')),varargin));
varargin(Idx:Idx+1) = [];
% Compile list of multi-Gaussian models
multiGaussModels = cell(maxGaussians,1);
multiGaussModels{1} = @dd_onegauss;
if maxGaussians>=2, multiGaussModels{2} = @dd_twogauss;   end
if maxGaussians>=3, multiGaussModels{3} = @dd_threegauss; end
if maxGaussians>=4, multiGaussModels{4} = @dd_fourgauss;  end
if maxGaussians>=5, multiGaussModels{5} = @dd_fivegauss;  end
for i = 6:maxGaussians
    multiGaussModels{i} =  mixmodels({multiGaussModels{i-1},@dd_onegauss});
end

%If the user has specified som boundaries then set the models boundaries appropiately
if ~isempty(Upper) || ~isempty(Lower)
    %Run over all multigauss models
    for i=1:maxGaussians
        %Get the info about the models
        info = multiGaussModels{i}();
        boundary = zeros(1,info.nparam);
        ParamNames = {info.parameters(:).name};
        %Find the parameter indexes for the FWHM and mean distances
        FWHMidx = contains(ParamNames,'FWHM');
        Distidx = contains(ParamNames,'distance');
        if ~isempty(Upper)
            boundary(FWHMidx) = Upper(2); %FWHM upper bound
            boundary(Distidx) = Upper(1); %Mean distaces upper bound
            boundary(~(Distidx | FWHMidx)) = 1; %Amplitudes upper bound
            if numel(Upper)>2
                boundary(numel(boundary)+1) = Upper(3); %Modulation depth upper bound
                boundary(numel(boundary)+1:numel(boundary) + numel(Upper(4:end))) = Upper(4:end); %Background upper bound
            end
            UpperBounds{i} = boundary;
        else
            UpperBounds = [];
        end
        boundary = zeros(1,info.nparam);
        if ~isempty(Lower)
            boundary(FWHMidx) = Lower(2); %FWHM lower bound
            boundary(Distidx) = Lower(1); %Mean distaces lower bound
            boundary(~(Distidx | FWHMidx)) = 0; %Amplitudes lower bound
            if numel(Lower)>2
                boundary(numel(boundary)+1) = Lower(3); %Modulation depth lower bound
                boundary(numel(boundary)+1:numel(boundary) + numel(Lower(4:end))) = Lower(4:end); %Background upper bound
            end
            LowerBounds{i} = boundary;
        else
            LowerBounds = [];
        end
    end
else
    %Otherwise just pass them empty to use the model defaults
    LowerBounds = [];
    UpperBounds = [];
end


if ~isempty(BckgModel)
    if isempty(LowerBounds)
        for i=1:maxGaussians
            info = multiGaussModels{i}();
            range = [info.parameters(:).range];
            Plower = range(1:2:end-1);
            infoB = BckgModel();
            range = [infoB.parameters(:).range];
            Blower = range(1:2:end-1);
            LowerBounds{i} = [Plower 0 Blower];
            
        end
    end
    if isempty(UpperBounds)
        for i=1:maxGaussians
            info = multiGaussModels{i}();
            range = [info.parameters(:).range];
            Pupper = range(2:2:end);
            infoB = BckgModel();
            range = [infoB.parameters(:).range];
            Bupper = range(2:2:end);
            UpperBounds{i} = [Pupper 0 Bupper];
        end
    end
    for i=1:maxGaussians
        DistModel = multiGaussModels{i};
        info = DistModel();
        Nparam = info.nparam;
        Pparam = [info.parameters(:).default];
        infoB = BckgModel();
        Bparam = [infoB.parameters(:).default];
        range = [infoB.parameters(:).range];
        Blower = range(1:2:end-1);
        Bupper = range(2:2:end);
        timeMultiGaussModels{i} = @(t,param) (1-param(Nparam+1) + param(Nparam+1)*dipolarkernel(t,r)*DistModel(r,param(1:Nparam)) ).*BckgModel(t,param(Nparam+2:end));
        param0{i} = [Pparam 0.5 Bparam];
    end
end

% Run fitting and model selection to see which multi-Gauss model is optimal
if ~isempty(BckgModel)
    [nGaussOpt,metrics,fitparams] = selectmodel(timeMultiGaussModels,S,t,method,param0,'Lower',LowerBounds,'Upper',UpperBounds,varargin);
else
    [nGaussOpt,metrics,fitparams] = selectmodel(multiGaussModels,S,r,K,method,'Lower',LowerBounds,'Upper',UpperBounds,varargin);
end

% Calculate the distance distribution for the optimal multi-Gauss model
param = fitparams{nGaussOpt};
optModel = multiGaussModels{nGaussOpt};
info = optModel();
Pfit = optModel(r,param(1:info.nparam));

if nargout>4
    Peval = zeros(maxGaussians,numel(r));
    for i=1:maxGaussians
        info = multiGaussModels{i}();
        p = fitparams{i};
        Peval(i,:) = multiGaussModels{i}(r,p(1:info.nparam));
    end
end

return
