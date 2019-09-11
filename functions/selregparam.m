%
% SELREGPARAM Selection of optimal regularization parameter
%
%   alpha = SELREGPARAM(S,K,'type','method')
%   Returns the optimal regularization parameter (alpha) from a range of
%   regularization parameter candidates (alphas). The parameter for the
%   regularization type given by ('type') is computed based on the input
%   signal (S), and the dipolar kernel (K).
%   The method employed for the selection of the regularization parameter
%   can be specified as the ('method') input argument.
%
%   alpha = SELREGPARAM(S,K,'type',{'method1',...,'methodN'})
%   If multiple selection methods are passed as a cell array of strings,
%   the function returns (alpha) as an N-point array of optimal
%   regularization parameters corresponding to the input methods,
%
%   [alpha,F,alphas] = SELREGPARAM(...,{'method1',...,'methodN'})
%   If requested, the second output argument returns the model selection
%   functionals of correspnding to the different input selection methods. A
%   third output argument (alphas) returns a vector with the alpha candidate
%   values evaluated in the search.
%
%   alpha = SELREGPARAM({S1,S2,...},{K1,K2,...},r,'type','method')
%   Passing multiple signals/kernels enables selection of the regularization
%   parameter for global fitting of the regularization model to a
%   single distribution. The global fit weights are automatically computed
%   according to their contribution to ill-posedness.
%
%   alpha = SELREGPARAM(...,'Property',Values)
%   Additional (optional) arguments can be passed as property-value pairs.
%
% The properties to be passed as options can be set in any order.
%
%   'NonNegConstrained' - True/false to enforce non-negativity (default=true)
%
%   'RegOrder' - Order of the regularization operator L (default = 2).
%
%   'Refine' - True/false to enforce a second search aroud the optimal value
%              with a finer grid to achieve a better value of the optimum.
%
%   'HuberParameter' - Huber parameter used in the 'huber' model (default = 1.35).
%
%   'GlobalWeights' - Array of weighting coefficients for the individual signals in
%                     global fitting regularization.
%
%   'TolFun' - Optimizer function tolerance.
%
%   'NoiseLevel' - Array of noise levels of the input signals.
%
%   'Range' - Range of alpha-value candidates to evaluate
%

% This file is a part of DeerAnalysis. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function [OptRegParam,Functionals,RegParamRange] = selregparam(S,K,RegType,SelectionMethod,varargin)

%--------------------------------------------------------------------------
% Parse & Validate Required Input
%--------------------------------------------------------------------------


%Turn off warnings to avoid ill-conditioned warnings 
warning('off','MATLAB:nearlySingularMatrix')

%Check if user requested some options via name-value input
[TolFun,NonNegConstrained,NoiseLevel,Refine,GlobalWeights,HuberParameter,RegParamRange,RegOrder] ...
    = parseoptional({'TolFun','NonNegConstrained','NoiseLevel','Refine','GlobalWeights','HuberParameter','Range','RegOrder'},varargin);

if ~iscell(S)
    S = {S};
end
if ~iscell(K)
    K = {K};
end
if isempty(RegOrder)
    RegOrder = 2;
else
    validateattributes(RegOrder,{'numeric'},{'scalar','nonnegative'})
end
if length(K)~=length(S)
    error('The number of kernels and signals must be equal.')
end
if ~isempty(GlobalWeights)
    validateattributes(GlobalWeights,{'numeric'},{'nonnegative'})
    if length(GlobalWeights) ~= length(S)
        error('The same number of global fit weights as signals must be passed.')
    end
    if sum(GlobalWeights)~=1
        error('The sum of the global fit weights must equal 1.')
    end
end
for i=1:length(S)
    if ~iscolumn(S{i})
        S{i} = S{i}.';
    end
    if ~isreal(S{i})
        S{i} = real(S{i});
    end
    if length(S{i})~=size(K{i},1)
        error('K and signal arguments must fulfill size(K,1)==length(S).')
    end
    validateattributes(S{i},{'numeric'},{'nonempty'},mfilename,'S')
    validateattributes(K{i},{'numeric'},{'nonempty'},mfilename,'K')
end
%Validate the selection methods input
allowedMethodInputs = {'lr','lc','dp','cv','gcv','rgcv','srgcv','aic','bic','aicc','rm','ee','ncp','gml','mcl'};
if iscell(SelectionMethod)
    for i=1:length(SelectionMethod)
        if strcmp(SelectionMethod{i},'all')
            SelectionMethod = allowedMethodInputs;
            break;
        end
        validateattributes(SelectionMethod{i},{'char'},{'nonempty'})
        SelectionMethod{i} = validatestring(SelectionMethod{i},allowedMethodInputs);
    end
else
    validateattributes(SelectionMethod,{'char'},{'nonempty'})
    if strcmp(SelectionMethod,'all')
        SelectionMethod = allowedMethodInputs;
    else
        SelectionMethod = validatestring(SelectionMethod,allowedMethodInputs);
        SelectionMethod = {SelectionMethod};
    end
end

%--------------------------------------------------------------------------
% Parse & Validate Optional Input
%--------------------------------------------------------------------------
warning('off','all')

%Validate nonNegLSQsolTol input
if isempty(TolFun)
    TolFun = 1e-9;
else
    validateattributes(TolFun,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'nonNegLSQsolTol')
end
%Validate RegType input
if isempty(RegType)
    RegType = 'tikhonov';
else
    validateattributes(RegType,{'char'},{'nonempty'})
    allowedInput = {'tikhonov','tv','huber','hubercirc','berhu'};
    RegType = validatestring(RegType,allowedInput);
end
if isempty(HuberParameter)
    HuberParameter = 1.35;
else
    validateattributes(HuberParameter,{'numeric'},{'scalar','nonempty','nonnegative'})
end
%Validate NoiseLevel input
if isempty(NoiseLevel)
    for i=1:length(S)
        NoiseLevel(i) = noiselevel(S{i});
    end
else
    if length(NoiseLevel)~=length(S)
        error('The same number of noise levels as signals must be provided.')
    end
    validateattributes(NoiseLevel,{'numeric'},{'scalar','nonempty','nonnegative'})
end
%Validate NonNegConstrained input
if isempty(NonNegConstrained)
    NonNegConstrained = true;
else
    validateattributes(NonNegConstrained,{'logical'},{'nonempty'},mfilename,'NonNegConstrained')
end


%--------------------------------------------------------------------------
% Preparations
%--------------------------------------------------------------------------
nr = size(K{1},2);
%Get regularization operator
L = regoperator(nr,RegOrder);
%Get range of potential alpha values candidates
if isempty(RegParamRange)
    RegParamRange = regparamrange(K{1},L);
else
    validateattributes(RegParamRange,{'numeric'},{'nonempty','nonnegative'},mfilename,'RegParamRange')
end
%Update number of points just to make sure
nPoints = length(RegParamRange);
%Initialize arrays
Residual = zeros(1,nPoints);
Penalty = zeros(1,nPoints);
HuberParameterSet = zeros(1,nPoints);
%Initialize cells
PseudoInverse = cell(1,nPoints);
P = cell(1,nPoints);
InfluenceMatrix = cell(1,nPoints);

%Preallocate common matrix multiplication
% KtK = K.'*K;


%--------------------------------------------------------------------------
% Pseudo-Inverses and Ps
%--------------------------------------------------------------------------

for i=1:nPoints %Loop over all regularization parameter values
    
    [Q,KtS,weights] = lsqcomponents(S,K,L,RegParamRange(i),RegType);
    InitialGuess = zeros(nr,1);
    if NonNegConstrained
        P{i} = fnnls(Q,KtS,InitialGuess,TolFun);
    else
        P{i}  = Q\KtS;
    end
    for idx=1:length(S)
        Q = lsqcomponents(S{idx},K{idx},L,RegParamRange(i),'tikhonov');
        PseudoInverse{idx,i} = Q\K{idx}.';
        switch lower(RegType)
            case 'tikhonov'
                Penalty(idx,i) = 1/sqrt(2)*norm(L*P{i});
            case 'tv'
                Penalty(idx,i) = sum(sqrt((L*P{i}).^2 + 1e-24));
            case 'huber'
                Penalty(idx,i) = sum(sqrt((L*P{i}/HuberParameter).^2 + 1 ) - 1);
        end
        Residual(idx,i) = 1/sqrt(2)*norm(K{idx}*P{i} - S{idx});
        InfluenceMatrix{idx,i} = K{idx}*PseudoInverse{idx,i};
    end
end

%In case variables are in matrix form reshape them to vectors
Functionals = cell(length(SelectionMethod),1);
OptRegParam = zeros(length(SelectionMethod),1);
OptHuberParam = zeros(length(SelectionMethod),1);

%--------------------------------------------------------------------------
% Selection methods for optimal regularization parameter
%--------------------------------------------------------------------------

%If multiple selection methods are requested then process them sequentially
for MethodIndex = 1:length(SelectionMethod)
    Functional = zeros(1,nPoints);
    for SIndex = 1:length(S)
        nr = length(S{SIndex});
        switch lower(SelectionMethod{MethodIndex})
            
            case 'lr' %L-curve Minimum-Radius method (LR)
                Eta = log(Penalty(SIndex,:));
                Rho = log(Residual(SIndex,:));
                Functional = Functional + weights(SIndex)*((((Rho - min(Rho))/(max(Rho) - min(Rho))).^2 + ((Eta - min(Eta))/(max(Eta) - min(Eta))).^2));
                
            case 'lc' %L-curve Maximum-Curvature method (LC)
                d1Residual = gradient(log(Residual(SIndex,:)));
                d2Residual = gradient(d1Residual);
                d1Penalty = gradient(log(Penalty(SIndex,:)));
                d2Penalty = gradient(d1Penalty);
                Functional = Functional + weights(SIndex)*((d1Residual.*d2Penalty - d2Residual.*d1Penalty)./(d1Residual.^2 + d1Penalty.^2).^(3/2));
                
            case 'dp' %Discrepancy principle (DP)
                SafetyFactor = 1;
                Index = Residual(SIndex,:)/sqrt(nr) <= SafetyFactor*NoiseLevel(SIndex);
                Functional(Index) = Functional(Index) - weights(SIndex)*(RegParamRange(Index));
                
            case 'cv' %Cross validation (CV)
                for i=1:nPoints
                    InfluenceDiagonal = diag(InfluenceMatrix{SIndex,i});
                    Functional(i) = Functional(i) + weights(SIndex)*(sum(abs(S{SIndex} - K{SIndex}*(P{i})./(ones(nr,1) - InfluenceDiagonal)).^2));
                end
                
            case 'gcv' %Generalized Cross Validation (GCV)
                for i=1:nPoints
                    Functional(i) = Functional(i) + weights(SIndex)*(Residual(i)^2/((1 - trace(InfluenceMatrix{SIndex,i})/nr)^2));
                end
                
            case 'rgcv' %Robust Generalized Cross Validation (rGCV)
                TuningParameter = 0.9;
                for i=1:nPoints
                    Functional(i) = Functional(i) + weights(SIndex)*(Residual(SIndex,i)^2/((1 - trace(InfluenceMatrix{SIndex,i})/nr)^2)*(TuningParameter + (1 - TuningParameter)*trace(InfluenceMatrix{SIndex,i}^2)/nr));
                end
                
            case 'srgcv' %Strong Robust Generalized Cross Validation (srGCV)
                TuningParameter = 0.8;
                for i=1:nPoints
                    Functional(i) = Functional(i) + weights(SIndex)*(Residual(SIndex,i)^2/((1 - trace(InfluenceMatrix{SIndex,i})/nr)^2)*(TuningParameter + (1 - TuningParameter)*trace(PseudoInverse{SIndex,i}'*PseudoInverse{SIndex,i})/nr));
                end
                
            case 'aic' %Akaike information criterion (AIC)
                Criterion = 2;
                for i=1:nPoints
                    Functional(i) = Functional(i) + weights(SIndex)*(nr*log(Residual(SIndex,i)^2/nr) + Criterion*trace(InfluenceMatrix{SIndex,i}));
                end
                
            case 'bic' %Bayesian information criterion (BIC)
                Criterion = log(nr);
                for i=1:nPoints
                    Functional(i) = Functional(i) + weights(SIndex)*(nr*log(Residual(SIndex,i)^2/nr) + Criterion*trace(InfluenceMatrix{SIndex,i}));
                end
                
            case 'aicc' %Corrected Akaike information criterion (AICC)
                for i=1:nPoints
                    Criterion = 2*nr/(nr-trace(InfluenceMatrix{SIndex,i})-1);
                    Functional(i) = Functional(i) + weights(SIndex)*(nr*log(Residual(SIndex,i)^2/nr) + Criterion*trace(InfluenceMatrix{SIndex,i}));
                end
                
            case 'rm' %Residual method (RM)
                for i=1:nPoints
                    Scaling = K{SIndex}.'*(eye(size(InfluenceMatrix{SIndex,i})) - InfluenceMatrix{SIndex,i});
                    Functional(i) = Functional(i) + weights(SIndex)*(Residual(SIndex,i)^2/sqrt(trace(Scaling'*Scaling)));
                end
                
            case 'ee' %Extrapolated Error (EE)
                for i=1:nPoints
                    Functional(i) = Functional(i) + weights(SIndex)*(Residual(SIndex,i)^2/norm(K{SIndex}.'*(K{SIndex}*P{i} - S{SIndex})));
                end
                
            case 'ncp' %Normalized Cumulative Periodogram (NCP)
                for i=1:nPoints
                    ResidualPeriodogram = abs(fft(K{SIndex}*P{i} - S{SIndex})).^2;
                    WhiteNoisePowerSpectrum = zeros(length(ResidualPeriodogram),1);
                    ResidualPowerSpectrum = zeros(length(ResidualPeriodogram),1);
                    for j=1:length(ResidualPeriodogram) - 1
                        ResidualPowerSpectrum(j)  = norm(ResidualPeriodogram(2:j+1),1)/norm(ResidualPeriodogram(2:end),1);
                        WhiteNoisePowerSpectrum(j) = j/(length(ResidualPeriodogram) - 1);
                    end
                    Functional(i) = Functional(i) + weights(SIndex)*(norm(ResidualPowerSpectrum - WhiteNoisePowerSpectrum));
                end
                
            case 'gml' %Generalized Maximum Likelihood (GML)
                Treshold = 1e-9;
                for i=1:nPoints
                    try %Once crushed beacause of eig(NaN)
                        EigenValues = eig(eye(size(InfluenceMatrix{SIndex,i})) - InfluenceMatrix{SIndex,i});
                    catch
                        EigenValues = 0;
                    end
                    EigenValues(EigenValues < Treshold) = 0;
                    NonZeroEigenvalues = real(EigenValues(EigenValues~=0));
                    Functional(i) = Functional(i) + weights(SIndex)*(S{SIndex}'*(S{SIndex} - K{SIndex}*P{i})/nthroot(prod(NonZeroEigenvalues),length(NonZeroEigenvalues)));
                end
                
            case 'mcl' %Mallows' C_L (MCL)
                for i=1:nPoints
                    Functional(i) = Residual(i)^2 + 2*NoiseLevel(SIndex)^2*trace(InfluenceMatrix{SIndex,i}) - 2*nr*NoiseLevel(SIndex)^2;
                end
                
        end
    end
    
    %Get optimal index of the selection functionals
    [~,Index] = min(Functional);
    %Store the optimal regularization parameter
    OptRegParam(MethodIndex) = RegParamRange(Index);
    OptHuberParam(MethodIndex) = HuberParameterSet(Index);
    Functionals{MethodIndex} =  Functional;
    
end

%If requested refine the grid search in a more precise search
if Refine
    RefineLength = 10;
    varargin{end+1} = 'Refine';
    varargin{end+1} =  false;
    OptIndex = 10;
    while any(OptIndex == RefineLength) || any(OptIndex == 1)
        if OptIndex == RefineLength
            FineRegParamRange = linspace(0.5*OptRegParam(1),2*OptRegParam(1),RefineLength);
        elseif OptIndex == 1
            FineRegParamRange = linspace(OptRegParam(1),2*OptRegParam(1),RefineLength);
        else
            FineRegParamRange = linspace(0.5*OptRegParam(1),2*OptRegParam(1),RefineLength);
        end
        varargin{end+1} = 'Range';
        varargin{end+1} = FineRegParamRange;
        [RefinedOptRegParam,RefinedFunctionals] = selregparam(S,K,RegType,SelectionMethod,varargin);
        for i=1:length(Functionals)
            Functionals{i} = [Functionals{i} RefinedFunctionals{i}];
        end
        RegParamRange = [RegParamRange FineRegParamRange];
        OptIndex = find(FineRegParamRange == RefinedOptRegParam(1));
        OptRegParam  = RefinedOptRegParam;
    end
end

%Turn warnings back on
warning('on','MATLAB:nearlySingularMatrix')

end

