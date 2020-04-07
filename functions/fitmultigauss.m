%
% FITMULTIGAUSS Multi-Gauss fitting of a distance distribution
%
%   P = FITMULTIGAUSS(V,K,r,Nmax,method)
%   Fits a multi-Gauss parametric distance distribution model to the dipolar
%   signal (V), using the dipolar kernel (K) and distance axis (r). The function
%   compares multi-Gaussian distributions with up to a maximum number of Gaussians
%   given by (Nmax) and determines the optimum one using the model selection
%   criterion given in (method) ('AIC', 'BIC', or 'AICc'). The fitted
%   distribution is returned in P.
%
%   P = FITMULTIGAUSS(V,t,r,Nmax,method)
%   If a the default kernel is to be used, the time axis (t) can be passed
%   instead of the kernel.
%
%   P = FITMULTIGAUSS({V1,V1,...},{K1,K2,K3},r,Nmax,method)
%   Passing multiple signals/kernels enables distance-domain global fitting
%   of the parametric models to single distributions. 
%   The multiple signals are passed as a cell array of arrays of sizes N1,N2,...
%   and a cell array of kernel matrices with sizes N1xM,N2xM,... must be 
%   passed as well.
%
%   P = FITMULTIGAUSS({V1,V1,...},{t1,t2,t3},r,Nmax,method)
%   Similarly, time-domain global fitting can be used when passing time-domain
%   and the model time axes {t1,t2,...} of the corresponding signals.
%
%   [P,param,Pci,paramci,opt,metrics,Peval] = FITMULTIGAUSS(___)
%   If requested alongside the distribution (P), the optimal fit model
%   parameters (param), the optimal number of Gaussians (opt) and
%   evaluated selection metrics (metrics) are returned. The fitted distance
%   distributions fitted for ech multigauss model can be requested as a
%   fifth output argument (Peval).
%
%   P = FITMULTIGAUSS(...,'Property',Value)
%   Additional (optional) arguments can be passed as property-value pairs.
%
%   'Background' - Function handle to the background model to be fitted
%                 along the multigauss distance distribution model.
%                 Requires the time-axis to be passed instead of the
%                 kernel. Can be used with global fitting, where the same
%                 model will be applied to all signals.
%
%   'Lower' - Array [<r>_min FWHM_min] containing the lower bound for the
%             FWHM and mean distance of all the Gaussians.
%   'Upper' -  Array [<r>_max FWHM_max] containing the upper bound values
%              for the FWHM and mean distance of all the Gaussians.
%
%   See "help fitparamodel" for a detailed list of other property-value pairs
%   accepted by the function.
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function [Pfit,param,Pfitci,paramci,nGaussOpt,metrics,Peval] = fitmultigauss(Vs,Ks,r,maxGaussians,method,varargin)


% Validate user input (S, K, r, and method are validated in lower-level functions)
if nargin<4
    error('Not enough input arguments.')
else
    validateattributes(maxGaussians,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'maxGaussians')
end
if nargin<5
    method = 'aicc';
elseif nargin==6
    varargin = [{method} varargin];
    method = 'aicc';
end

% Parse the optional parameters in the varargin
optionalProperties = {'Upper','Lower','Background','internal::parselater'};
[Upper,Lower,BckgModel] = parseoptional(optionalProperties,varargin);

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

% Remove used options from varargin so they are not passed to fitparamodel
for i = 1:numel(optionalProperties)
    Idx = find(cellfun(@(x)(ischar(x) && strcmpi(x,optionalProperties{i})),varargin));
    varargin(Idx:Idx+1) = [];
end

if nargin>3
    % Include internal option in order for fitparamodel to return the covariance matrix
    varargin = [varargin {'internal::returncovariancematrix'} {true}];
end



%Parse the required inputs for global fitting
if ~iscell(Vs)
   Vs = {Vs}; 
end
if ~iscell(Ks)
   Ks = {Ks}; 
end
for i=1:numel(Ks)
    if ~all(size(Ks{i}) > 1)
        ts{i} = Ks{i};
        Ks{i} = dipolarkernel(ts{i},r);
    end
end
if ~isempty(BckgModel) && ~exist('ts','var')
    error('Time axes must be provided for a time-domain fit.')
end
Nsignals = numel(Vs);

% Compile list of multi-Gaussian models
multiGaussModels = cell(maxGaussians,1);
multiGaussModels{1} = @dd_gauss;
if maxGaussians>=2, multiGaussModels{2} = @dd_gauss2; end
if maxGaussians>=3, multiGaussModels{3} = @dd_gauss3; end
if maxGaussians>=4, multiGaussModels{4} = @dd_gauss4; end
if maxGaussians>=5, multiGaussModels{5} = @dd_gauss5; end
for i = 6:maxGaussians
    multiGaussModels{i} =  mixmodels(multiGaussModels{i-1},@dd_gauss);
end

% If the user has specified som boundaries then set the models boundaries appropiately
if ~isempty(Upper) || ~isempty(Lower)
    % Run over all multi-Gauss models
    for i = 1:maxGaussians
        % Get the info about the models
        info = multiGaussModels{i}();
        boundary = zeros(1,info.nparam);
        ParamNames = {info.parameters(:).name};
        % Find the parameter indexes for the FWHM and mean distances
        FWHMidx = contains(ParamNames,'FWHM');
        Distidx = contains(ParamNames,'distance');
        if ~isempty(Upper)
            boundary(FWHMidx) = Upper(2); % FWHM upper bound
            boundary(Distidx) = Upper(1); % Mean distaces upper bound
            boundary(~(Distidx | FWHMidx)) = 1; % Amplitudes upper bound
            if numel(Upper)>2
                boundary(numel(boundary)+1) = Upper(3); % Modulation depth upper bound
                boundary(numel(boundary)+1:numel(boundary) + numel(Upper(4:end))) = Upper(4:end); % Background upper bound
            end
            UpperBounds{i} = boundary;
        else
            UpperBounds = [];
        end
        boundary = zeros(1,info.nparam);
        if ~isempty(Lower)
            boundary(FWHMidx) = Lower(2); % FWHM lower bound
            boundary(Distidx) = Lower(1); % Mean distaces lower bound
            boundary(~(Distidx | FWHMidx)) = 0; % Amplitudes lower bound
            if numel(Lower)>2
                boundary(numel(boundary)+1) = Lower(3); % Modulation depth lower bound
                boundary(numel(boundary)+1:numel(boundary) + numel(Lower(4:end))) = Lower(4:end); %Background upper bound
            end
            LowerBounds{i} = boundary;
        else
            LowerBounds = [];
        end
    end
else
    % Otherwise just pass them empty to use the model defaults
    LowerBounds = [];
    UpperBounds = [];
end


if ~isempty(BckgModel)
    if isempty(LowerBounds)
        for i = 1:maxGaussians
            info = multiGaussModels{i}();
            range = [info.parameters(:).range];
            Plower = range(1:2:end-1);
            infoB = BckgModel();
            range = [infoB.parameters(:).range];
            Blower = range(1:2:end-1);
            LowerBounds{i} = [Plower repmat([0 Blower],1,Nsignals)];
            
        end
    end
    if isempty(UpperBounds)
        for i = 1:maxGaussians
            info = multiGaussModels{i}();
            range = [info.parameters(:).range];
            Pupper = range(2:2:end);
            infoB = BckgModel();
            range = [infoB.parameters(:).range];
            Bupper = range(2:2:end);
            UpperBounds{i} = [Pupper repmat([1 Bupper],1,Nsignals)];
        end
    end
    for i = 1:maxGaussians
        DistModel = multiGaussModels{i};
        info = DistModel();
        Nparam = info.nparam;
        Pparam = [info.parameters(:).default];
        infoB = BckgModel();
        Bparam = [infoB.parameters(:).default];
        lampars = Nparam + (1+numel(Bparam))*(1:Nsignals)-numel(Bparam);
        Bpars = lampars + 1;
        lam0 = 0.25;
        timeMultiGaussModels{i} = @(t,param,idx) (1 - param(lampars(idx)) + param(lampars(idx))*dipolarkernel(t,r)*DistModel(r,param(1:Nparam)) ).*BckgModel(t,param(Bpars(idx):Bpars(idx)+numel(Bparam)-1));
        param0{i} = [Pparam repmat([lam0 Bparam],1,Nsignals)];
    end
end

% Run fitting and model selection to see which multi-Gauss model is optimal
if ~isempty(BckgModel)
    [nGaussOpt,metrics,fitparams,paramcis] = selectmodel(timeMultiGaussModels,Vs,ts,method,param0,'Lower',LowerBounds,'Upper',UpperBounds,varargin);
else
    [nGaussOpt,metrics,fitparams,paramcis] = selectmodel(multiGaussModels,Vs,r,Ks,method,'Lower',LowerBounds,'Upper',UpperBounds,varargin);
end

% Calculate the distance distribution for the optimal multi-Gauss model
param = fitparams{nGaussOpt};
paramci = paramcis{nGaussOpt}{1};
optModel = multiGaussModels{nGaussOpt};
info = optModel();
Pfit = optModel(r,param(1:info.nparam));

if nargin>3
    % Compute Jacobian for current multi-Gauss model
    covmatrix = paramcis{nGaussOpt}{2};
    critical = paramcis{nGaussOpt}{3};
    covmatrix = covmatrix(1:info.nparam,1:info.nparam);
    jacobian = jacobianest(@(par)optModel(r,par),param(1:info.nparam));
    % Calculate the confidence bands for the distance distribution
    modelvariance = arrayfun(@(idx)full(jacobian(idx,:))*covmatrix*full(jacobian(idx,:)).',1:numel(r)).';
    upper = Pfit + critical*sqrt(modelvariance);
    lower = Pfit - critical*sqrt(modelvariance);
    upper = max(upper,0);
    lower = max(lower,0);
    Pfitci = [upper(:) lower(:)];
end

if nargout>6
    Peval = zeros(maxGaussians,numel(r));
    for i = 1:maxGaussians
        info = multiGaussModels{i}();
        p = fitparams{i};
        Peval(i,:) = multiGaussModels{i}(r,p(1:info.nparam));
    end
end

return
