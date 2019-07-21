function [OptRegParam,Functionals,Lcurve,OptHuberParam] = selectRegParam(RegParamRange,Signal,Kernel,RegMatrix,SelectionMethod,varargin)

%--------------------------------------------------------------------------
% Parse & Validate Required Input
%--------------------------------------------------------------------------
if ~iscolumn(Signal)
    Signal = Signal';
end
validateattributes(Signal,{'numeric'},{'nonempty'},mfilename,'Signal')
validateattributes(RegParamRange,{'numeric'},{'nonempty','nonnegative'},mfilename,'RegParamRange')
validateattributes(Kernel,{'numeric'},{'nonempty'},mfilename,'Kernel')
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
        validatestring(SelectionMethod{i},allowedMethodInputs);
    end
else
    validateattributes(SelectionMethod,{'char'},{'nonempty'})
    if strcmp(SelectionMethod,'all')
        SelectionMethod = allowedMethodInputs;
    else
        validatestring(SelectionMethod,allowedMethodInputs);
        SelectionMethod = {SelectionMethod};
    end
end

%--------------------------------------------------------------------------
% Parse & Validate Optional Input
%--------------------------------------------------------------------------
%Check if user requested some options via name-value input
[nonNegLSQsolTol,NonNegConstrained,RegType,NoiseLevel] = parseOptional({'nonNegLSQsolTol','NonNegConstrained','RegType','NoiseLevel'},varargin);

warning('off','all')

%Validate nonNegLSQsolTol input
if isempty(nonNegLSQsolTol)
    nonNegLSQsolTol = 1e-9;
else
    validateattributes(nonNegLSQsolTol,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'nonNegLSQsolTol')
end
%Validate RegType input
if isempty(RegType)
    RegType = 'tikhonov';
else
    validateattributes(RegType,{'char'},{'nonempty'})
    allowedInput = {'tikhonov','tv','huber','hubercirc','berhu'};
    validatestring(RegType,allowedInput);
end
%Validate NoiseLevel input
if isempty(NoiseLevel)
    NoiseLevel = std(Signal(end-floor(0.25*length(Signal)):end));
else
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
DipolarDimension = length(Signal);

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
Percentile = 80:2.5:95;

%--------------------------------------------------------------------------
% Pseudo-Inverses and Distributions
%--------------------------------------------------------------------------

for i=1:nPoints %Loop over all regularization parameter values
    
    switch RegType
        
        %--------------------------------------------------------------------------
        % Total variation (L1) Penalty
        %--------------------------------------------------------------------------
        case 'tv'
            
            %Unconstrained distributions required for construction of correct pseudoinverse
            Distribution{i} = zeros(DipolarDimension,1);
            %Define options for fsolve
            Solveroptions = optimoptions(@fsolve,'Display','off','Algorithm','trust-region-reflective');
            Distribution{i} = fsolve(@(Distribution) Kernel'*(Kernel*Distribution - Signal) + ...
                RegParamRange(i)^2*(RegMatrix)'*((RegMatrix*Distribution)./sqrt((RegMatrix*Distribution).^2 + 1e-24)) ...
                ,Distribution{i},Solveroptions);
            
            %Compute pseudoinverse by analytical expression (only defined for unconstrained distribution)
            localTVPseudoInverse = RegMatrix'*((RegMatrix./sqrt((RegMatrix*Distribution{i}).^2 + 1e-24)));
            PseudoInverse{i} = (Kernel'*Kernel + RegParamRange(i)^2*localTVPseudoInverse)\Kernel';
            
            %If constrained solution is requested
            if NonNegConstrained
                %Now that we have the pseudoinverse, the constrained distribution can be computed quickly
                Distribution{i} = fnnls(Kernel'*Kernel + RegParamRange(i)^2*localTVPseudoInverse,Kernel'*Signal,[],nonNegLSQsolTol);
            end
            
            %Get last variables required for the selection functionals
            Penalty(i) = norm(RegMatrix*Distribution{i});
            %Residual defined this way so that later Residual(i)^2 equals 1/2*norm(...)^2
            Residual(i) = 1/sqrt(2)*norm(Kernel*Distribution{i} - Signal);
            InfluenceMatrix{i} = Kernel*PseudoInverse{i};
            
            %--------------------------------------------------------------------------
            % Pseudo-Huber (L2-L1) Penalty
            %--------------------------------------------------------------------------
        case 'huber'
            
            %Compute approximate optimal Huber Parameter for current value of regularization parameter
            if i>1
                    %(Skip the first value since it requires previous unconstrained distributions)
                    Solveroptions = optimoptions(@fsolve,'Display','off','Algorithm','trust-region-reflective');
                    HuberParameter = fsolve(@(HuberParameter) Kernel'*(Kernel*PseudoInverse{i-1}*Signal - Signal) + ...
                        RegParamRange(i-1)^2/HuberParameter^2*(RegMatrix')*((RegMatrix*PseudoInverse{i-1}*Signal)./sqrt((RegMatrix*PseudoInverse{i-1}*Signal/HuberParameter).^2 + 1)) ...
                        ,1.35,Solveroptions);
                    % Ensure that Huber parameter does not get negative (can lead to crashes later on)
                    HuberParameter = abs(HuberParameter);
                    HuberParameterSet(i) = HuberParameter;
            else
               HuberParameter = 1.35; 
            end
            %Unconstrained distributions required for construction of correct pseudoinverse
            Distribution{i} = zeros(DipolarDimension,1);
            try
                %Define options for fsolve
                Solveroptions = optimoptions(@fsolve,'Display','off','Algorithm','trust-region-reflective');
                Distribution{i} = fsolve(@(Distribution) Kernel'*(Kernel*Distribution - Signal) + ...
                    RegParamRange(i)^2/(HuberParameter^2)*(RegMatrix')*((RegMatrix*Distribution)./sqrt((RegMatrix*(Distribution)/HuberParameter).^2 + 1)) ...
                    ,Distribution{i},Solveroptions);
                
                % Compute pseudoinverse by analytical expression (only defined for unconstrained distribution)
                HuberTerm = 1/(HuberParameter^2)*((RegMatrix)'*(RegMatrix./sqrt((RegMatrix*Distribution{i}/HuberParameter).^2 + 1)));
                PseudoInverse{i} = (Kernel'*Kernel + RegParamRange(i)^2*HuberTerm)\Kernel';
                
                % If constrained solution is requested
                if NonNegConstrained
                    %Now that we have the pseudoinverse, the constrained distribution can be computed quickly
                    Distribution{i} = fnnls(Kernel'*Kernel + RegParamRange(i)^2*HuberTerm,Kernel'*Signal);
                end
            catch
                % Worst-case scenario when fsolve crashes, use this bad method to get at least something
                fprintf('Huber Selection: Failed \n')
                PseudoInverseTemp = (Kernel'*Kernel + RegParamRange(i)^2*(RegMatrix')*RegMatrix)\Kernel';
                DistributionTemp  = PseudoInverseTemp*Signal;
                HuberTerm = ((RegMatrix/HuberParameter)'*diag(1./sqrt((RegMatrix*DistributionTemp/HuberParameter).^2 + 1))*(RegMatrix/HuberParameter));
                PseudoInverse{i} = (Kernel'*Kernel + RegParamRange(i)^2*HuberTerm)\Kernel';
                Distribution{i} = fnnls(Kernel'*Kernel + RegParamRange(i)^2*HuberTerm,Kernel'*Signal);
            end
            
            %Get last variables required for the selection functionals
            Penalty(i) = sum(sqrt((RegMatrix*Distribution{i}/HuberParameter).^2 + 1 ) - 1);
            Residual(i) = 1/sqrt(2)*norm(Kernel*Distribution{i} - Signal);
            InfluenceMatrix{i} = Kernel*(PseudoInverse{i});
            
            %--------------------------------------------------------------------------
            % Circular-Symmetric Huber Penalty
            %--------------------------------------------------------------------------
        case 'hubercirc'
            PseudoInverseTemp = (Kernel'*Kernel + RegParamRange(i)^2*(RegMatrix')*RegMatrix)\Kernel';
            if options.unconstrainedProblem
                DistributionTemp  = PseudoInverseTemp*Signal; %unconstrained Tikhonov solution
            else
                DistributionTemp = fnnls(Kernel'*Kernel + RegParamRange(i)^2*(RegMatrix')*RegMatrix,Kernel'*Signal ); %non-negative Tikhonov solution
            end
            HuberParameter = norm(RegMatrix*DistributionTemp);
            HuberParameter = HuberParameter + [eps ];
            for j=1:length(HuberParameter)
                
                %First case of circ. sym. Huber function ( Norm < Huber )
                if norm(RegMatrix*DistributionTemp) < HuberParameter(j)
                    PseudoInverse{j,i} = (Kernel'*Kernel + RegParamRange(i)^2*(RegMatrix')*RegMatrix)\Kernel';
                    if NonNegConstrained
                        Distribution{j,i} = fnnls(Kernel'*Kernel + RegParamRange(i)^2*(RegMatrix')*RegMatrix,Kernel'*Signal);
                    else
                        Distribution{j,i}  = PseudoInverse{j,i}*Signal;
                    end
                    
                else
                    
                    %Second case of circ. sym. Huber function ( Norm > Huber )
                    PseudoInverse{j,i} = (Kernel'*Kernel + RegParamRange(i)^2*2*HuberParameter(j)/norm(RegMatrix*DistributionTemp)*(RegMatrix*DistributionTemp)'*RegMatrix)\Kernel';
                    
                    if i>1
                        try
                            Solveroptions = optimoptions(@fsolve,'Display','off','Algorithm','trust-region-reflective');
                            HuberParameter(j) = fsolve(@(HuberParameter) Kernel'*(Kernel*Distribution{j,i-1} - Signal) + ...
                                RegParamRange(i)^2*2*HuberParameter/norm(RegMatrix*Distribution{j,i-1})*(RegMatrix*Distribution{j,i-1})'*RegMatrix ,mean(edgeIndicator(RegMatrix*DistributionTemp,'gradient')),Solveroptions);
                            % Ensure that Huber parameter does not get negative (can lead to crushed later on)
                            HuberParameter(j) = abs(HuberParameter(j));
                        catch
                        end
                    end
                    Distribution{j,i}  =  PseudoInverse{j,i}*Signal;
                    %Define options for independent-fsolve
                    Solveroptions = optimoptions(@fsolve,'Display','off','Algorithm','trust-region-reflective');
                    Distribution{j,i} = fsolve(@(Distribution) Kernel'*(Kernel*Distribution - Signal) + ...
                        RegParamRange(i)^2*2*HuberParameter(j)/norm(RegMatrix*Distribution)*(RegMatrix*Distribution)'*RegMatrix,Distribution{j,i},Solveroptions);
                    
                    % Compute pseudoinverse by analytical expression (only defined for unconstrained distribution)
                    HuberTerm = 2*HuberParameter(j)/norm(RegMatrix*Distribution{j,i})*(RegMatrix*Distribution{j,i})'*RegMatrix;
                    PseudoInverse{j,i} = (Kernel'*Kernel + RegParamRange(i)^2*HuberTerm)\Kernel';
                    
                    % If constrained solution is requested
                    if ~options.unconstrainedProblem
                        %Now that we have the pseudoinverse, the constrained distribution can be computed quickly
                        Distribution{i} = fnnls(Kernel'*Kernel + RegParamRange(i)^2*HuberTerm,Kernel'*Signal);
                    end
                    
                    PseudoInverse{j,i} = (Kernel'*Kernel + RegParamRange(i)^2*2*HuberParameter(j)/norm(RegMatrix*Distribution{j,i})*(RegMatrix')*RegMatrix)\Kernel';
                end
                
                Penalty(j,i) = norm(RegMatrix*Distribution{j,i});
                Residual(j,i) = 1/sqrt(2)*norm(Kernel*Distribution{j,i} - Signal);
                InfluenceMatrix{j,i} = Kernel*PseudoInverse{j,i};
            end
            
            %--------------------------------------------------------------------------
            % Inverse Huber Berhu (L1-L2) Penalty
            %--------------------------------------------------------------------------
        case 'berhu'
            HuberParameter = options.HuberParameter;
            for j=1:length(HuberParameter)
                PseudoInverse{j,i} = (Kernel'*Kernel + RegParamRange(i)^2*(RegMatrix')*RegMatrix)\Kernel';
                Distribution{j,i}  = PseudoInverse{j,i}*Signal; %unconstrained Tikhonov solution
                HuberParameter = prctile(edgeIndicator(DistributionTemp,RegMatrix,'fire'),Percentile);
                Auxiliary = RegMatrix*Distribution{j,i};
                Temp = (RegMatrix)'*diag(1./sqrt(1e-24 + (Auxiliary.^2)))*(RegMatrix);
                BerhuGradient = zeros(length(Auxiliary)+options.Derivative,1);
                for k=1:length(Auxiliary)
                    if abs(Auxiliary(k)) > HuberParameter(j)
                        BerhuGradient(k) = Auxiliary(k)/HuberParameter(j);
                    else
                        BerhuGradient(k) =  Temp(k);
                    end
                end
                PseudoInverse{j,i} = (Kernel'*Kernel + RegParamRange(i)^2*BerhuGradient)\Kernel';
                if options.unconstrainedProblem
                    Distribution{j,i}  = PseudoInverse{j,i}*Signal;
                else
                    Distribution{j,i} = fnnls(Kernel'*Kernel + RegParamRange(i)^2*BerhuGradient,Kernel'*Signal);
                end
                Penalty(j,i) = norm(Kernel*Distribution{j,i} - Signal);
                Residual(j,i) = 1/sqrt(2)*norm(Kernel*Distribution{j,i} - Signal);
                InfluenceMatrix{j,i} = Kernel*PseudoInverse{j,i};
            end
            
            %--------------------------------------------------------------------------
            % Tikhonov (L2) Penalty
            %--------------------------------------------------------------------------
        case 'tikhonov'
            PseudoInverse{i} = (Kernel'*Kernel + RegParamRange(i)^2*(RegMatrix')*RegMatrix)\Kernel';
            if NonNegConstrained
                Distribution{i} = fnnls(Kernel'*Kernel + RegParamRange(i)^2*(RegMatrix')*RegMatrix,Kernel'*Signal,[],nonNegLSQsolTol); %non-negative Tikhonov solution
            else
                Distribution{i}  = PseudoInverse{i}*Signal; %unconstrained Tikhonov solution
            end
            Penalty(i) = 1/sqrt(2)*norm(RegMatrix*Distribution{i});
            Residual(i) = 1/sqrt(2)*norm(Kernel*Distribution{i} - Signal);
            InfluenceMatrix{i} = Kernel*PseudoInverse{i};
    end
end

%In case variables are in matrix form reshape them to vectors
Distribution = reshape(Distribution,1,[]);
PseudoInverse = reshape(PseudoInverse,1,[]);
Residual = reshape(Residual,1,[]);
InfluenceMatrix = reshape(InfluenceMatrix,1,[]);
Penalty = reshape(Penalty,1,[]);
nPoints = length(Distribution);
Functionals = cell(length(SelectionMethod));
OptRegParam = zeros(length(SelectionMethod),1);
OptHuberParam = zeros(length(SelectionMethod),1);

%--------------------------------------------------------------------------
% Selection methods for optimal regularization parameter
%--------------------------------------------------------------------------

%If multiple selection methods are requested then process them sequentially
for MethodIndex = 1:length(SelectionMethod)
    
    switch SelectionMethod{MethodIndex}
        
        case 'lr' %L-curve Minimum-Radius method (LR)
            Eta = log(Penalty);
            Rho = log(Residual);
            Functional = ((Rho - min(Rho))/(max(Rho) - min(Rho))).^2 + ((Eta - min(Eta))/(max(Eta) - min(Eta))).^2;
            
        case 'lc' %L-curve Maximum-Curvature method (LC)
            d1Residual = gradient(log(Residual));
            d2Residual = gradient(d1Residual);
            d1Penalty = gradient(log(Penalty));
            d2Penalty = gradient(d1Penalty);
            Functional = (d1Residual.*d2Penalty - d2Residual.*d1Penalty)./(d1Residual.^2 + d1Penalty.^2).^(3/2);
            
        case 'dp' %Discrepancy principle (DP)
            SafetyFactor = 1;
            Index = Residual/sqrt(DipolarDimension) <= SafetyFactor*NoiseLevel;
            Functional = -RegParamRange(Index);
            
        case 'cv' %Cross validation (CV)
            for i=1:nPoints
                InfluenceDiagonal = diag(InfluenceMatrix{i});
                Functional(i) = sum(abs(Signal - Kernel*(Distribution{i})./(ones(DipolarDimension,1) - InfluenceDiagonal)).^2);
            end
            
        case 'gcv' %Generalized Cross Validation (GCV)
            for i=1:nPoints
                Functional(i) = Residual(i)^2/((1 - trace(InfluenceMatrix{i})/DipolarDimension)^2);
            end
            
        case 'rgcv' %Robust Generalized Cross Validation (rGCV)
            TuningParameter = 0.9;
            for i=1:nPoints
                Functional(i) = Residual(i)^2/((1 - trace(InfluenceMatrix{i})/DipolarDimension)^2)*(TuningParameter + (1 - TuningParameter)*trace(InfluenceMatrix{i}^2)/DipolarDimension);
            end
            
        case 'srgcv' %Strong Robust Generalized Cross Validation (srGCV)
            TuningParameter = 0.8;
            for i=1:nPoints
                Functional(i) = Residual(i)^2/((1 - trace(InfluenceMatrix{i})/DipolarDimension)^2)*(TuningParameter + (1 - TuningParameter)*trace(PseudoInverse{i}'*PseudoInverse{i})/DipolarDimension);
            end
            
        case 'aic' %Akaike information criterion (AIC)
            Criterion = 2;
            for i=1:nPoints
                Functional(i) = DipolarDimension*log(Residual(i)^2/DipolarDimension) + Criterion*trace(InfluenceMatrix{i});
            end
            
        case 'bic' %Bayesian information criterion (BIC)
            Criterion = log(DipolarDimension);
            for i=1:nPoints
                Functional(i) = DipolarDimension*log(Residual(i)^2/DipolarDimension) + Criterion*trace(InfluenceMatrix{i});
            end
            
        case 'aicc' %Corrected Akaike information criterion (AICC)
            for i=1:nPoints
                Criterion = 2*DipolarDimension/(DipolarDimension-trace(InfluenceMatrix{i})-1);
                Functional(i) = DipolarDimension*log(Residual(i)^2/DipolarDimension) + Criterion*trace(InfluenceMatrix{i});
            end
            
        case 'rm' %Residual method (RM)
            for i=1:nPoints
                Scaling = Kernel'*(eye(size(InfluenceMatrix{i})) - InfluenceMatrix{i});
                Functional(i) = Residual(i)^2/sqrt(trace(Scaling'*Scaling));
            end
            
        case 'ee' %Extrapolated Error (EE)
            for i=1:nPoints
                Functional(i) = Residual(i)^2/norm(Kernel'*(Kernel*Distribution{i} - Signal));
            end
            
        case 'ncp' %Normalized Cumulative Periodogram (NCP)
            for i=1:nPoints
                ResidualPeriodogram = abs(fft(Kernel*Distribution{i} - Signal)).^2;
                WhiteNoisePowerSpectrum = zeros(length(ResidualPeriodogram),1);
                ResidualPowerSpectrum = zeros(length(ResidualPeriodogram),1);
                for j=1:length(ResidualPeriodogram) - 1
                    ResidualPowerSpectrum(j)  = norm(ResidualPeriodogram(2:j+1),1)/norm(ResidualPeriodogram(2:end),1);
                    WhiteNoisePowerSpectrum(j) = j/(length(ResidualPeriodogram) - 1);
                end
                Functional(i) = norm(ResidualPowerSpectrum - WhiteNoisePowerSpectrum);
            end
            
        case 'gml' %Generalized Maximum Likelihood (GML)
            Treshold = 1e-9;
            for i=1:nPoints
                try %Once crushed beacause of eig(NaN)
                    EigenValues = eig(eye(size(InfluenceMatrix{i})) - InfluenceMatrix{i});
                catch
                    EigenValues = 0;
                end
                EigenValues(EigenValues < Treshold) = 0;
                NonZeroEigenvalues = real(EigenValues(EigenValues~=0));
                Functional(i) = Signal'*(Signal - Kernel*Distribution{i})/nthroot(prod(NonZeroEigenvalues),length(NonZeroEigenvalues));
            end
            
        case 'mcl' %Mallows' C_L (MCL)
            for i=1:nPoints
                Functional(i) = Residual(i)^2 + 2*NoiseLevel^2*trace(InfluenceMatrix{i}) - 2*DipolarDimension*NoiseLevel^2;
            end
            
    end
    
    %Get optimal index of the selection functionals
    [~,Index] = min(Functional);
    %Store the optimal regularization parameter
    OptRegParam(MethodIndex) = RegParamRange(Index);
    OptHuberParam(MethodIndex) = HuberParameterSet(Index);
    if nargout>1
        Functionals{MethodIndex} =  Functional;
    end
    
end

%If requested get the L-curve parameters and return them
if nargout>2
    Eta = log(Penalty);
    Rho = log(Residual);
    Lcurve = [Eta Rho];
end


end

