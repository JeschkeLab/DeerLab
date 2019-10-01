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
%
%   'Lower' - Array [<r>_min FWHM_min] containing the lower bound for the 
%   FWHM and mean distance of all the Gaussians.
%
%   'Upper' -  Array [<r>_max FWHM_max] containing the upper bound values 
%   for the FWHM and mean distance of all the Gaussians.
%
%   See "help fitparamodel" for a detailed list of other property-value pairs
%   accepted by the function.
%

% This file is a part of DeerAnalysis. License is MIT (see LICENSE.md).
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

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

warning('off','DA:parseoptional')
%Parse the optional parameters in the varargin
[Upper,Lower] = parseoptional({'Upper','Lower'},varargin);
warning('on','DA:parseoptional')
if isempty(Upper)
elseif ~isempty(Upper) && length(Upper)~=2
    error('Upper property must be an array [<r>_max FWHM_max]')
end
if ~isempty(Lower) && length(Lower)~=2
    error('Lower property must be an array [<r>_min FWHM_min]')
end
%Remove the Lower and Upper options from varargin so they are not passed to fitparamodel
Idx = find(cellfun(@(x)(ischar(x) && contains(x,'Upper')),varargin));
varargin(Idx:Idx+1) = [];
Idx = find(cellfun(@(x)(ischar(x) && contains(x,'Lower')),varargin));
varargin(Idx:Idx+1) = [];

% Compile list of multi-Gaussian models
multiGaussModels = cell(maxGaussians,1);
multiGaussModels{1} = @rd_onegaussian;
if maxGaussians>=2, multiGaussModels{2} = @rd_twogaussian;   end
if maxGaussians>=3, multiGaussModels{3} = @rd_threegaussian; end
if maxGaussians>=4, multiGaussModels{4} = @rd_fourgaussian;  end
if maxGaussians>=5, multiGaussModels{5} = @rd_fivegaussian;  end
for i = 6:maxGaussians
    multiGaussModels{i} =  mixmodels({multiGaussModels{i-1},@rd_onegaussian});
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
            UpperBounds{i} = boundary;
        else
            UpperBounds = [];
        end
        if ~isempty(Lower)
            boundary(FWHMidx) = Lower(2); %FWHM lower bound
            boundary(Distidx) = Lower(1); %Mean distaces lower bound
            boundary(~(Distidx | FWHMidx)) = 0; %Amplitudes lower bound
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


% Run fitting and model selection to see which multi-Gauss model is optimal
[nGaussOpt,metrics,fitparams] = selectmodel(multiGaussModels,S,r,K,method,'Lower',LowerBounds,'Upper',UpperBounds,varargin);

% Calculate the distance distribution for the optimal multi-Gauss model
param = fitparams{nGaussOpt};
optModel = multiGaussModels{nGaussOpt};
Pfit = optModel(r,param);

if nargout>4
    Peval = zeros(maxGaussians,numel(r));
    for i=1:maxGaussians
        Peval(i,:) = multiGaussModels{i}(r,fitparams{i});
    end
end

return
