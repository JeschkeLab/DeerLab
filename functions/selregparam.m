%
%  SELREGPARAM Selection of optimal regularization parameter
%
%    alpha = SELREGPARAM(S,K,r,'type','method')
%    Returns the optimal regularization parameter (alpha) from a range of
%    regularization parameter candidates (alphas). The parameter for the
%    regularization type given by ('type') is computed based on the input
%    signal (S), and the dipolar kernel (K) on a distance axis (r).
%    The method employed for the selection of the regularization parameter
%    can be specified as the ('method') input argument.
%
%    alpha = SELREGPARAM(S,K,r,'type',{'method1',...,'methodN'})
%    If multiple selection methods are passed as a cell array of strings,
%    the function returns (alpha) as an N-point array of optimal
%    regularization parameters corresponding to the input methods,
%
%    [alpha,F,alphas] = SELREGPARAM(...,{'method1',...,'methodN'})
%    If requested, the second output argument returns the model selection
%    functionals of correspnding to the different input selection methods. A
%    third output argument (alphas) returns a vector with the alpha candidate
%    values evaluated in the search.
%
%    alpha = SELREGPARAM({S1,S2,...},{K1,K2,...},r,'type','method')
%    Passing multiple signals/kernels enables selection of the regularization
%    parameter for global fitting of the regularization model to a
%    single distribution. The global fit weights are automatically computed
%    according to their contribution to ill-posedness.
%
%    alpha = SELREGPARAM(...,'Property',Values)
%    Additional (optional) arguments can be passed as property-value pairs.
%
%  The properties to be passed as options can be set in any order.
%
%    'NonNegConstrained' - True/false to enforce non-negativity (default=true)
%
%    'RegOrder' - Order of the regularization operator L (default = 2).
%
%    'Refine' - True/false to enforce a second search aroud the optimal value
%               with a finer grid to achieve a better value of the optimum.
%
%    'HuberParameter' - Huber parameter used in the 'huber' model (default = 1.35).
%
%    'GlobalWeights' - Array of weighting coefficients for the individual signals in
%                      global fitting regularization.
%
%    'TolFun' - Optimizer function tolerance.
%
%    'NoiseLevel' - Array of noise levels of the input signals.
%
%    'Range' - Range of alpha-value candidates to evaluate
%

%  This file is a part of DeerAnalysis. License is MIT (see LICENSE.md).
%  Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function [alphaOpt,Functionals,alphaRange,Residuals,Penalties] = selregparam(S,K,r,RegType,SelectionMethod,varargin)

%--------------------------------------------------------------------------
%  Parse & Validate Required Input
%--------------------------------------------------------------------------


% Turn off warnings to avoid ill-conditioned warnings
warning('off','MATLAB:nearlySingularMatrix')

% Check if user requested some options via name-value input
[TolFun,NonNegConstrained,NoiseLevel,GlobalWeights,HuberParameter,alphaRange,RegOrder,Search] ...
    = parseoptional({'TolFun','NonNegConstrained','NoiseLevel','GlobalWeights','HuberParameter','Range','RegOrder','Search'},varargin);

if ~iscell(S)
    S = {S};
end
if ~iscell(K)
    K = {K};
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
        S{i} = error('Input signal(s) cannot be complex.');
    end
    if length(S{i})~=size(K{i},1)
        error('K and signal arguments must fulfill size(K,1)==length(S).')
    end
    validateattributes(S{i},{'numeric'},{'nonempty'},mfilename,'S')
    validateattributes(K{i},{'numeric'},{'nonempty'},mfilename,'K')
end
% Validate the selection methods input
allowedMethodInputs = {'lr','lc','cv','gcv','rgcv','srgcv','aic','bic','aicc','rm','ee','ncp','gml','mcl'};
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
%  Parse & Validate Optional Input
%--------------------------------------------------------------------------
warning('off','all')

if isempty(RegOrder)
    RegOrder = 2;
else
    validateattributes(RegOrder,{'numeric'},{'scalar','nonnegative'})
end

% Validate nonNegLSQsolTol input
if isempty(TolFun)
    TolFun = 1e-9;
else
    validateattributes(TolFun,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'nonNegLSQsolTol')
end
% Validate RegType input
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
% Validate NoiseLevel input
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
% Validate NonNegConstrained input
if isempty(NonNegConstrained)
    NonNegConstrained = true;
else
    validateattributes(NonNegConstrained,{'logical'},{'nonempty'},mfilename,'NonNegConstrained')
end
% Validate Search input
if isempty(Search)
    Search = 'golden';
else
    validateattributes(Search,{'char'},{'nonempty'},mfilename,'Search')
    Search = validatestring(Search,{'golden','exhaustive'});
end

%--------------------------------------------------------------------------
%  Preparations
%--------------------------------------------------------------------------
nr = size(K{1},2);
% Get regularization operator
L = regoperator(nr,RegOrder);
% Get range of potential alpha values candidates
if isempty(alphaRange)
    alphaRange = regparamrange(K{1},L);
else
    validateattributes(alphaRange,{'numeric'},{'nonempty','nonnegative'},mfilename,'RegParamRange')
end
% Update number of points just to make sure
nPoints = length(alphaRange);
% Initialize arrays
HuberParameterSet = zeros(1,nPoints);



%--------------------------------------------------------------------------
%  Golden search algorithm
%--------------------------------------------------------------------------
switch  lower(Search)
    
    case 'golden'
        
        for jj = 1:numel(SelectionMethod)
            intervalStart = min(log(alphaRange));
            intervalEnd = max(log(alphaRange));
            %         intervalStart = -10;
            %         intervalEnd = 10;
            epsilon = 0.1;               
            iter = 10;                      
            tau = double((sqrt(5)-1)/2);    
            k = 0;                          
            logalpha1 = intervalStart + (1-tau)*(intervalEnd-intervalStart);      
            logalpha2 = intervalStart + tau*(intervalEnd-intervalStart);
            
            fcnval1 = evalalpha(exp(logalpha1),SelectionMethod(jj));
            fcnval2 = evalalpha(exp(logalpha2),SelectionMethod(jj));
            
            alphasEvaluated = [logalpha1 logalpha2];
            Functional = [fcnval1 fcnval2];
            
            while abs(intervalEnd-intervalStart) > epsilon && k < iter
                k = k + 1;
                if(fcnval1<fcnval2)
                    intervalEnd = logalpha2;
                    logalpha2 = logalpha1;
                    logalpha1 = intervalStart + (1 - tau)*(intervalEnd - intervalStart);
                    
                    fcnval1 = evalalpha(exp(logalpha1),SelectionMethod(jj));
                    fcnval2 = evalalpha(exp(logalpha2),SelectionMethod(jj));
                    
                    alphasEvaluated(end+1) = logalpha1;
                    Functional(end+1) = fcnval1;
                else
                    intervalStart = logalpha1;
                    logalpha1 = logalpha2;
                    logalpha2 = intervalStart+tau*(intervalEnd - intervalStart);
                    
                    fcnval1 = evalalpha(exp(logalpha1),SelectionMethod(jj));
                    fcnval2 = evalalpha(exp(logalpha2),SelectionMethod(jj));
                    
                    alphasEvaluated(end+1) = logalpha2;
                    Functional(end+1) = fcnval2;
                end
                k=k+1;
            end
            Functionals{jj} = Functional;
            if(fcnval1<fcnval2)
                alphaOpt(jj) = exp(logalpha1);
            else
                alphaOpt(jj) = exp(logalpha2);
            end
            alphaRanges{jj} = alphasEvaluated;
            
        end
        a = 1;
        %--------------------------------------------------------------------------
        %  Exhaustive search
        %--------------------------------------------------------------------------
    case 'exhaustive'
        
        for ii=1:numel(alphaRange)
            Functional(ii,:) = evalalpha(alphaRange(ii),SelectionMethod);
        end
        for jj = 1:numel(SelectionMethod)        
        % Get optimal index of the selection functionals
        [~,Index] = min(Functional(:,jj));
        % Store the optimal regularization parameter
        alphaOpt(jj) = alphaRange(Index);
        Functionals{jj} =  Functional(:,jj);
        end
end

% Turn warnings back on
warning('on','MATLAB:nearlySingularMatrix')


%--------------------------------------------------------------------------
    function Functional = evalalpha(alpha,SelectionMethod)
        
        
        %--------------------------------------------------------------------------
        %  Pseudo-Inverses and Ps
        %--------------------------------------------------------------------------
        
        [Q,KtS,weights] = lsqcomponents(S,r,K,L,alpha,RegType,HuberParameter,GlobalWeights);
        InitialGuess = zeros(numel(r),1);
        if NonNegConstrained
            P = fnnls(Q,KtS,InitialGuess,TolFun);
        else
            P  = Q\KtS;
        end
        for idx = 1:length(S)
            Q = lsqcomponents(S{idx},r,K{idx},L,alpha,'tikhonov',HuberParameter,GlobalWeights);
            PseudoInverse{idx} = Q\K{idx}.';
            switch lower(RegType)
                case 'tikhonov'
                    Penalty(idx) = norm(L*P);
                case 'tv'
                    Penalty(idx) = sum(sqrt((L*P).^2 + 1e-24));
                case 'huber'
                    Penalty(idx) = sum(sqrt((L*P/HuberParameter).^2 + 1 ) - 1);
            end
            Residual(idx) = norm(K{idx}*P - S{idx});
            InfluenceMatrix{idx} = K{idx}*PseudoInverse{idx};
        end
        
        %--------------------------------------------------------------------------
        %  Selection methods for optimal regularization parameter
        %--------------------------------------------------------------------------
        
        % If multiple selection methods are requested then process them sequentially
            Functional = zeros(length(SelectionMethod),1);
        for MethodIndex = 1:length(SelectionMethod)
            for SIndex = 1:length(S)
                nr = length(S{SIndex});
                switch lower(SelectionMethod{MethodIndex})
                    
                    case 'lr' % L-curve Minimum-Radius method (LR)
                        Eta = log(Penalty(SIndex));
                        Rho = log(Residual(SIndex));
                        Functional(MethodIndex) = Functional(MethodIndex) + weights(SIndex)*((((Rho - min(Rho))/(max(Rho) - min(Rho))).^2 + ((Eta - min(Eta))/(max(Eta) - min(Eta))).^2));
                        
                    case 'lc' % L-curve Maximum-Curvature method (LC)
                        d1Residual = gradient(log(Residual(SIndex)));
                        d2Residual = gradient(d1Residual);
                        d1Penalty = gradient(log(Penalty(SIndex)));
                        d2Penalty = gradient(d1Penalty);
                        Functional(MethodIndex) = Functional(MethodIndex) + weights(SIndex)*((d1Residual.*d2Penalty - d2Residual.*d1Penalty)./(d1Residual.^2 + d1Penalty.^2).^(3/2));
                        
                    case 'cv' % Cross validation (CV)
                        InfluenceDiagonal = diag(InfluenceMatrix{SIndex});
                        Functional(MethodIndex) = Functional(MethodIndex) + weights(SIndex)*(sum(abs(S{SIndex} - K{SIndex}*(P)./(ones(nr,1) - InfluenceDiagonal)).^2));
                        
                    case 'gcv' % Generalized Cross Validation (GCV)
                        Functional(MethodIndex) = Functional(MethodIndex) + weights(SIndex)*(Residual(SIndex)^2/((1 - trace(InfluenceMatrix{SIndex})/nr)^2));
                        
                    case 'rgcv' % Robust Generalized Cross Validation (rGCV)
                        TuningParameter = 0.9;
                        Functional(MethodIndex) = Functional(MethodIndex) + weights(SIndex)*(Residual(SIndex)^2/((1 - trace(InfluenceMatrix{SIndex})/nr)^2)*(TuningParameter + (1 - TuningParameter)*trace(InfluenceMatrix{SIndex}^2)/nr));
                        
                    case 'srgcv' % Strong Robust Generalized Cross Validation (srGCV)
                        TuningParameter = 0.8;
                        Functional(MethodIndex) = Functional(MethodIndex) + weights(SIndex)*(Residual(SIndex)^2/((1 - trace(InfluenceMatrix{SIndex})/nr)^2)*(TuningParameter + (1 - TuningParameter)*trace(PseudoInverse{SIndex}'*PseudoInverse{SIndex})/nr));
                        
                    case 'aic' % Akaike information criterion (AIC)
                        Criterion = 2;
                        Functional(MethodIndex) = Functional(MethodIndex) + weights(SIndex)*(nr*log(Residual(SIndex)^2/nr) + Criterion*trace(InfluenceMatrix{SIndex}));
                        
                    case 'bic' % Bayesian information criterion (BIC)
                        Criterion = log(nr);
                        Functional(MethodIndex) = Functional(MethodIndex) + weights(SIndex)*(nr*log(Residual(SIndex)^2/nr) + Criterion*trace(InfluenceMatrix{SIndex}));
                        
                    case 'aicc' % Corrected Akaike information criterion (AICC)
                        Criterion = 2*nr/(nr-trace(InfluenceMatrix{SIndex})-1);
                        Functional(MethodIndex) = Functional(MethodIndex) + weights(SIndex)*(nr*log(Residual(SIndex)^2/nr) + Criterion*trace(InfluenceMatrix{SIndex}));
                        
                    case 'rm' % Residual method (RM)
                        Scaling = K{SIndex}.'*(eye(size(InfluenceMatrix{SIndex})) - InfluenceMatrix{SIndex});
                        Functional(MethodIndex) = Functional(MethodIndex) + weights(SIndex)*(Residual(SIndex)^2/sqrt(trace(Scaling'*Scaling)));
                        
                    case 'ee' % Extrapolated Error (EE)
                        Functional(MethodIndex) = Functional(MethodIndex) + weights(SIndex)*(Residual(SIndex)^2/norm(K{SIndex}.'*(K{SIndex}*P - S{SIndex})));
                        
                    case 'ncp' % Normalized Cumulative Periodogram (NCP)
                        ResidualPeriodogram = abs(fft(K{SIndex}*P - S{SIndex})).^2;
                        WhiteNoisePowerSpectrum = zeros(length(ResidualPeriodogram),1);
                        ResidualPowerSpectrum = zeros(length(ResidualPeriodogram),1);
                        for j=1:length(ResidualPeriodogram) - 1
                            ResidualPowerSpectrum(j)  = norm(ResidualPeriodogram(2:j+1),1)/norm(ResidualPeriodogram(2:end),1);
                            WhiteNoisePowerSpectrum(j) = j/(length(ResidualPeriodogram) - 1);
                        end
                        Functional(MethodIndex) = Functional(MethodIndex) + weights(SIndex)*(norm(ResidualPowerSpectrum - WhiteNoisePowerSpectrum));
                        
                    case 'gml' % Generalized Maximum Likelihood (GML)
                        Treshold = 1e-9;
                        EigenValues = eig(eye(size(InfluenceMatrix{SIndex})) - InfluenceMatrix{SIndex});
                        EigenValues(EigenValues < Treshold) = 0;
                        NonZeroEigenvalues = real(EigenValues(EigenValues~=0));
                        Functional(MethodIndex) = Functional(MethodIndex) + weights(SIndex)*(S{SIndex}'*(S{SIndex} - K{SIndex}*P)/nthroot(prod(NonZeroEigenvalues),length(NonZeroEigenvalues)));
                        
                    case 'mcl' % Mallows' C_L (MCL)
                        Functional(MethodIndex) = Functional(MethodIndex) +  weights(SIndex)*(Residual^2 + 2*NoiseLevel(SIndex)^2*trace(InfluenceMatrix{SIndex}) - 2*nr*NoiseLevel(SIndex)^2);
                        
                end
            end

        end
        
        
    end


end