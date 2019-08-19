%
% SELREGPARAM Selection of optimal regularization parameter
%
%   alpha = SELREGPARAM(alphas,S,K,L,'type','method')
%   Returns the optimal regularization parameter (alpha) from a range of
%   regularization parameter candidates (alphas). The parameter for the
%   regularization type given by ('type') is computed based on the input
%   signal (S), the dipolar kernel (K) and the regularization operator (L).
%   The method employed for the selection of the regularization parameter
%   can be specified as the ('method') input argument.
%
%   alpha = SELREGPARAM(alphas,S,K,L,'type',{'method1',...,'methodN'})
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
%   alpha = SELREGPARAM(alphas,{S1,S2,...},{K1,K2,...},r,L,'type','method')
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
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function [OptRegParam,Functionals,RegParamRange] = selregparam(RegParamRange,Signal,Kernel,RegMatrix,RegType,SelectionMethod,varargin)

%--------------------------------------------------------------------------
% Parse & Validate Required Input
%--------------------------------------------------------------------------

%Check if user requested some options via name-value input
[TolFun,NonNegConstrained,NoiseLevel,Refine,GlobalWeights,HuberParameter] = parseoptional({'TolFun','NonNegConstrained','NoiseLevel','Refine','GlobalWeights','HuberParameter'},varargin);

if ~iscell(Signal)
    Signal = {Signal};
end
if ~iscell(Kernel)
    Kernel = {Kernel};
end
if length(Kernel)~=length(Signal)
    error('The number of kernels and signals must be equal.')
end
if ~isempty(GlobalWeights)
    validateattributes(GlobalWeights,{'numeric'},{'nonnegative'})
    if length(GlobalWeights) ~= length(Signal)
        error('The same number of global fit weights as signals must be passed.')
    end
    if sum(GlobalWeights)~=1
        error('The sum of the global fit weights must equal 1.')
    end
end
for i=1:length(Signal)
    if ~iscolumn(Signal{i})
        Signal{i} = Signal{i}.';
    end
    if ~isreal(Signal{i})
        Signal{i} = real(Signal{i});
    end
    if length(Signal{i})~=size(Kernel{i},1)
        error('Kernel and signal arguments must fulfill size(Kernel,1)==length(Signal).')
    end
    validateattributes(Signal{i},{'numeric'},{'nonempty'},mfilename,'Signal')
    validateattributes(Kernel{i},{'numeric'},{'nonempty'},mfilename,'Kernel')
end
validateattributes(RegParamRange,{'numeric'},{'nonempty','nonnegative'},mfilename,'RegParamRange')
validateattributes(RegMatrix,{'numeric'},{'nonempty'},mfilename,'RegMatrix')
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
    for i=1:length(Signal)
        NoiseLevel(i) = noiselevel(Signal{i});
    end
else
    if length(NoiseLevel)~=length(Signal)
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
%Memoization
%--------------------------------------------------------------------------

persistent cachedData
if isempty(cachedData)
    cachedData =  java.util.LinkedHashMap;
end

hashKey = datahash({RegParamRange,Signal,Kernel,RegMatrix,SelectionMethod,varargin});
if cachedData.containsKey(hashKey)
    Output = cachedData.get(hashKey);
    [OptRegParam,Functionals,RegParamRange,OptHuberParam] = java2mat(Output);
    return
end

%--------------------------------------------------------------------------
% Preparations
%--------------------------------------------------------------------------
DipolarDimension = length(RegMatrix);

%Update number of points just to make sure
nPoints = length(RegParamRange);
%Initialize arrays
Residual = zeros(1,nPoints);
Penalty = zeros(1,nPoints);
Functional  = zeros(1,nPoints);
HuberParameterSet = zeros(1,nPoints);
%Initialize cells
PseudoInverse = cell(1,nPoints);
Distribution = cell(1,nPoints);
InfluenceMatrix = cell(1,nPoints);

%Preallocate common matrix multiplication
% KtK = Kernel.'*Kernel;


%--------------------------------------------------------------------------
% Pseudo-Inverses and Distributions
%--------------------------------------------------------------------------

for i=1:nPoints %Loop over all regularization parameter values
    
    [Q,KtS,weights] = lsqcomponents(Signal,Kernel,RegMatrix,RegParamRange(i),RegType);
    InitialGuess = zeros(DipolarDimension,1);
    if NonNegConstrained
        Distribution{i} = fnnls(Q,KtS,InitialGuess,TolFun);
    else
        Distribution{i}  = Q\KtS;
    end
    for idx=1:length(Signal)
        Q = lsqcomponents(Signal{idx},Kernel{idx},RegMatrix,RegParamRange(i),'tikhonov');
        PseudoInverse{idx,i} = Q\Kernel{idx}.';
        switch lower(RegType)
            case 'tikhonov'
                Penalty(idx,i) = 1/sqrt(2)*norm(RegMatrix*Distribution{i});
            case 'tv'
                Penalty(idx,i) = sum(sqrt((RegMatrix*Distribution{i}).^2 + 1e-24));
            case 'huber'
                Penalty(idx,i) = sum(sqrt((RegMatrix*Distribution{i}/HuberParameter).^2 + 1 ) - 1);
        end
        Residual(idx,i) = 1/sqrt(2)*norm(Kernel{idx}*Distribution{i} - Signal{idx});
        InfluenceMatrix{idx,i} = Kernel{idx}*PseudoInverse{idx,i};
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
    for SignalIndex = 1:length(Signal)
        DipolarDimension = length(Signal{SignalIndex});
        switch lower(SelectionMethod{MethodIndex})
            
            case 'lr' %L-curve Minimum-Radius method (LR)
                Eta = log(Penalty(SignalIndex,:));
                Rho = log(Residual(SignalIndex,:));
                Functional = Functional + weights(SignalIndex)*((((Rho - min(Rho))/(max(Rho) - min(Rho))).^2 + ((Eta - min(Eta))/(max(Eta) - min(Eta))).^2));
                
            case 'lc' %L-curve Maximum-Curvature method (LC)
                d1Residual = gradient(log(Residual(SignalIndex,:)));
                d2Residual = gradient(d1Residual);
                d1Penalty = gradient(log(Penalty(SignalIndex,:)));
                d2Penalty = gradient(d1Penalty);
                Functional = Functional + weights(SignalIndex)*((d1Residual.*d2Penalty - d2Residual.*d1Penalty)./(d1Residual.^2 + d1Penalty.^2).^(3/2));
                
            case 'dp' %Discrepancy principle (DP)
                SafetyFactor = 1;
                Index = Residual(SignalIndex,:)/sqrt(DipolarDimension) <= SafetyFactor*NoiseLevel(SignalIndex);
                Functional(Index) = Functional(Index) - weights(SignalIndex)*(RegParamRange(Index));
                
            case 'cv' %Cross validation (CV)
                for i=1:nPoints
                    InfluenceDiagonal = diag(InfluenceMatrix{SignalIndex,i});
                    Functional(i) = Functional(i) + weights(SignalIndex)*(sum(abs(Signal{SignalIndex} - Kernel{SignalIndex}*(Distribution{i})./(ones(DipolarDimension,1) - InfluenceDiagonal)).^2));
                end
                
            case 'gcv' %Generalized Cross Validation (GCV)
                for i=1:nPoints
                    Functional(i) = Functional(i) + weights(SignalIndex)*(Residual(i)^2/((1 - trace(InfluenceMatrix{SignalIndex,i})/DipolarDimension)^2));
                end
                
            case 'rgcv' %Robust Generalized Cross Validation (rGCV)
                TuningParameter = 0.9;
                for i=1:nPoints
                    Functional(i) = Functional(i) + weights(SignalIndex)*(Residual(SignalIndex,i)^2/((1 - trace(InfluenceMatrix{SignalIndex,i})/DipolarDimension)^2)*(TuningParameter + (1 - TuningParameter)*trace(InfluenceMatrix{SignalIndex,i}^2)/DipolarDimension));
                end
                
            case 'srgcv' %Strong Robust Generalized Cross Validation (srGCV)
                TuningParameter = 0.8;
                for i=1:nPoints
                    Functional(i) = Functional(i) + weights(SignalIndex)*(Residual(SignalIndex,i)^2/((1 - trace(InfluenceMatrix{SignalIndex,i})/DipolarDimension)^2)*(TuningParameter + (1 - TuningParameter)*trace(PseudoInverse{SignalIndex,i}'*PseudoInverse{SignalIndex,i})/DipolarDimension));
                end
                
            case 'aic' %Akaike information criterion (AIC)
                Criterion = 2;
                for i=1:nPoints
                    Functional(i) = Functional(i) + weights(SignalIndex)*(DipolarDimension*log(Residual(SignalIndex,i)^2/DipolarDimension) + Criterion*trace(InfluenceMatrix{SignalIndex,i}));
                end
                
            case 'bic' %Bayesian information criterion (BIC)
                Criterion = log(DipolarDimension);
                for i=1:nPoints
                    Functional(i) = Functional(i) + weights(SignalIndex)*(DipolarDimension*log(Residual(SignalIndex,i)^2/DipolarDimension) + Criterion*trace(InfluenceMatrix{SignalIndex,i}));
                end
                
            case 'aicc' %Corrected Akaike information criterion (AICC)
                for i=1:nPoints
                    Criterion = 2*DipolarDimension/(DipolarDimension-trace(InfluenceMatrix{SignalIndex,i})-1);
                    Functional(i) = Functional(i) + weights(SignalIndex)*(DipolarDimension*log(Residual(SignalIndex,i)^2/DipolarDimension) + Criterion*trace(InfluenceMatrix{SignalIndex,i}));
                end
                
            case 'rm' %Residual method (RM)
                for i=1:nPoints
                    Scaling = Kernel{SignalIndex}.'*(eye(size(InfluenceMatrix{SignalIndex,i})) - InfluenceMatrix{SignalIndex,i});
                    Functional(i) = Functional(i) + weights(SignalIndex)*(Residual(SignalIndex,i)^2/sqrt(trace(Scaling'*Scaling)));
                end
                
            case 'ee' %Extrapolated Error (EE)
                for i=1:nPoints
                    Functional(i) = Functional(i) + weights(SignalIndex)*(Residual(SignalIndex,i)^2/norm(Kernel{SignalIndex}.'*(Kernel{SignalIndex}*Distribution{i} - Signal{SignalIndex})));
                end
                
            case 'ncp' %Normalized Cumulative Periodogram (NCP)
                for i=1:nPoints
                    ResidualPeriodogram = abs(fft(Kernel{SignalIndex}*Distribution{i} - Signal{SignalIndex})).^2;
                    WhiteNoisePowerSpectrum = zeros(length(ResidualPeriodogram),1);
                    ResidualPowerSpectrum = zeros(length(ResidualPeriodogram),1);
                    for j=1:length(ResidualPeriodogram) - 1
                        ResidualPowerSpectrum(j)  = norm(ResidualPeriodogram(2:j+1),1)/norm(ResidualPeriodogram(2:end),1);
                        WhiteNoisePowerSpectrum(j) = j/(length(ResidualPeriodogram) - 1);
                    end
                    Functional(i) = Functional(i) + weights(SignalIndex)*(norm(ResidualPowerSpectrum - WhiteNoisePowerSpectrum));
                end
                
            case 'gml' %Generalized Maximum Likelihood (GML)
                Treshold = 1e-9;
                for i=1:nPoints
                    try %Once crushed beacause of eig(NaN)
                        EigenValues = eig(eye(size(InfluenceMatrix{SignalIndex,i})) - InfluenceMatrix{SignalIndex,i});
                    catch
                        EigenValues = 0;
                    end
                    EigenValues(EigenValues < Treshold) = 0;
                    NonZeroEigenvalues = real(EigenValues(EigenValues~=0));
                    Functional(i) = Functional(i) + weights(SignalIndex)*(Signal{SignalIndex}'*(Signal{SignalIndex} - Kernel{SignalIndex}*Distribution{i})/nthroot(prod(NonZeroEigenvalues),length(NonZeroEigenvalues)));
                end
                
            case 'mcl' %Mallows' C_L (MCL)
                for i=1:nPoints
                    Functional(i) = Residual(i)^2 + 2*NoiseLevel(SignalIndex)^2*trace(InfluenceMatrix{SignalIndex,i}) - 2*DipolarDimension*NoiseLevel(SignalIndex)^2;
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
        [RefinedOptRegParam,RefinedFunctionals] = selregparam(FineRegParamRange,Signal,Kernel,RegMatrix,RegType,SelectionMethod,varargin);
        for i=1:length(Functionals)
            Functionals{i} = [Functionals{i} RefinedFunctionals{i}];
        end
        RegParamRange = [RegParamRange FineRegParamRange];
        OptIndex = find(FineRegParamRange == RefinedOptRegParam(1));
        OptRegParam  = RefinedOptRegParam;
    end
end

%Store output result in the cache
Output = {OptRegParam,Functionals,RegParamRange,OptHuberParam};
cachedData = addcache(cachedData,hashKey,Output);


end

