function [Distribution,ConvergenceCurve] = obir(Signal,Kernel,DistanceAxis,RegType,RegMatrix,RegParam,varargin)
%--------------------------------------------------------------------------
% OSHER'S BREGMAN ITERATED REGULARIZATION (OBIR) METHOD
%--------------------------------------------------------------------------
% Regularization of DEER data by Bregman iterations [1] equiped with different
% penalty functionals of choice [2]. Bregman iterations allow the recovery of
% noiseless signal from the residuals of the previous iteration. The
% optimal solution is determined to be found when the residual standard
% deviation is equal or less than the noise standard deviation. Therefore
% the level of noise given as an input is determinant for the succesfull
% performance of the method. Nontheless, when the noise level can be well
% approximated, this method will yield better results than any other
% non-Bregman-iterated method even with optimal regularization parameter
% choice [3].
%
% Literature:
% [1] L.M. Bregman, USSR Computational Mathematics and Mathematical Physics (1967) 7, 3, 200-217
% [2] Charest et al, 2006 40th Annual Conference on Information Sciences and Systems
% [3] Yin et al, SIAM Journal on Imaging Sciences (2008) 1, 143-168
%--------------------------------------------------------------------------
% (C) 2019, Luis Fabregas, DeerAnalysis
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Parse & Validate Required Input
%--------------------------------------------------------------------------
if ~iscolumn(Signal)
    Signal = Signal';
end

%Get optional parameters
[NoiseLevelAim,Solver,MaxIter,TolFun,MaxFunEvals,DivergenceStop,MaxOuterIter,HuberParam,AxisHandle] = parseoptional({'NoiseLevelAim','Solver','MaxIter','TolFun','MaxFunEvals','DivergenceStop','MaxOuterIter','HuberParam','AxisHandle'},varargin);


if isempty(TolFun)
    TolFun = 1e-10;
else
    validateattributes(TolFun,{'numeric'},{'scalar','nonnegative'},mfilename,'TolFun')
end

if isempty(MaxOuterIter)
    MaxOuterIter = 5000;
else
    validateattributes(MaxOuterIter,{'numeric'},{'scalar','nonnegative'},mfilename,'MaxOuterIter')
end

if isempty(HuberParam)
    HuberParam = 1.35;
else
    validateattributes(HuberParam,{'numeric'},{'scalar','nonnegative'},mfilename,'HuberParam')
end

if isempty(MaxFunEvals)
    MaxFunEvals = 500000;
else
    validateattributes(MaxFunEvals,{'numeric'},{'scalar','nonnegative'},mfilename,'MaxFunEvals')
end

if isempty(MaxIter)
    MaxIter = 500000;
else
    validateattributes(MaxIter,{'numeric'},{'scalar','nonnegative'},mfilename,'MaxIter')
end

if isempty(Solver)
    Solver = 'fnnls';
else
    validateattributes(Solver,{'char'},{'nonempty'},mfilename,'Solver')
    validInputs = {'fnnls','fmincon'};
    validatestring(Solver,validInputs);
end

if isempty(NoiseLevelAim)
    NoiseLevelAim = noiselevel(Signal);
else
    validateattributes(NoiseLevelAim,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'NoiseLevelAim')
end
validateattributes(RegParam,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'RegParam')
validateattributes(Signal,{'numeric'},{'nonempty'},mfilename,'Signal')
validateattributes(Kernel,{'numeric'},{'nonempty'},mfilename,'Kernel')
validateattributes(RegMatrix,{'numeric'},{'nonempty'},mfilename,'RegMatrix')
checklengths(Signal,Kernel);

if nargin<3 || isempty(RegType)
    RegType = 'tikhonov';
elseif isa(RegType,'function_handle')
    RegFunctional = RegType;
    RegType = 'custom';
else
    validateattributes(RegType,{'char'},{'nonempty'},mfilename,'RegType')
    allowedInput = {'tikhonov','tv','huber'};
    validatestring(RegType,allowedInput);
end


%--------------------------------------------------------------------------
% Parse & Validate Optional Input
%--------------------------------------------------------------------------

if isempty(MaxOuterIter)
    MaxOuterIter = 200;
else
    validateattributes(MaxOuterIter,{'numeric'},{'scalar','nonempty'},mfilename,'MaxOuterIter')
end

if isempty(MaxIter)
    MaxIter = 200000;
else
    validateattributes(MaxIter,{'numeric'},{'scalar','nonempty'},mfilename,'MaxIter')
end


if isempty(MaxFunEvals)
    MaxFunEvals = 200000;
else
    validateattributes(MaxFunEvals,{'numeric'},{'scalar','nonempty'},mfilename,'MaxFunEvals')
end

if isempty(DivergenceStop)
    DivergenceStop = false;
else
    validateattributes(DivergenceStop,{'logical'},{'nonempty'},mfilename,'DivergenceStop')
end

%--------------------------------------------------------------------------
% Preparation
%--------------------------------------------------------------------------

%Initialize
SizeDistribution = length(Signal);
Subgradient = zeros(SizeDistribution,1);
Counter = 1;
Iteration = 1;
Distribution = zeros(SizeDistribution,1);

%--------------------------------------------------------------------------
% Osher's Bregman Iterations Algorithm
%--------------------------------------------------------------------------

Dimension = length(Signal);
InitialGuess = zeros(Dimension,1);
NonNegConst = zeros(Dimension,1);

while Iteration <= MaxOuterIter
    
    %Store privous iteration distribution
    CheckDistribution = Distribution;
    
    switch Solver
        case 'fmincon'
            %Define current minimization problem
            RegFunctional = regfunctional(RegType,Signal,RegMatrix,Kernel,RegParam,HuberParam);
            fminconFunctional = @(Distribution)OBIRFunctional(Distribution,RegFunctional,Subgradient);
            fminconOptions = optimset('GradObj','on','MaxFunEvals',MaxFunEvals,'Display','off','MaxIter',MaxIter);
            %Run minimzation
            Distribution =  fmincon(fminconFunctional,InitialGuess,[],[],[],[],NonNegConst,[],[],fminconOptions);
        case 'fnnls'
            
            %If using LSQ-based solvers then precompute the KtK and KtS input arguments
            KtS = Kernel.'*Signal - Subgradient;
            switch RegType
                case 'tikhonov'
                    %Constrained Tikhonov regularization
                    Q = (Kernel.'*Kernel) + RegParam^2*(RegMatrix.'*RegMatrix);
                case 'tv'
                    localDistribution = InitialGuess;
                    for j=1:500
                        prev = localDistribution;
                        %Compute pseudoinverse and unconst. distribution recursively
                        TVterm = RegMatrix'*((RegMatrix./sqrt((RegMatrix*localDistribution).^2 + 1e-24)));
                        localTVPseudoinv = (Kernel.'*Kernel + RegParam^2*TVterm)\Kernel.';
                        localDistribution = localTVPseudoinv*Signal;
                        change = norm(localDistribution - prev);
                        if round(change,5) ==0
                            break;
                        end
                    end
                    TVterm = RegMatrix.'*((RegMatrix./sqrt((RegMatrix*localDistribution).^2 + 1e-24)));
                    Q = (Kernel.'*Kernel + RegParam^2*TVterm);
                case 'huber'
                    localDistribution = InitialGuess;
                    for j=1:500
                        prev = localDistribution;
                        %Compute pseudoinverse and unconst. distribution recursively
                        HuberTerm = 1/(HuberParam^2)*((RegMatrix)'*(RegMatrix./sqrt((RegMatrix*localDistribution/HuberParam).^2 + 1)));
                        localHuberPseudoinv = (Kernel.'*Kernel + RegParam^2*HuberTerm)\Kernel.';
                        localDistribution = localHuberPseudoinv*Signal;
                        change = norm(localDistribution - prev);
                        if round(change,5) ==0
                            break;
                        end
                    end
                    HuberTerm = 1/(HuberParam^2)*((RegMatrix)'*(RegMatrix./sqrt((RegMatrix*localDistribution/HuberParam).^2 + 1)));
                    Q = (Kernel.'*Kernel + RegParam^2*HuberTerm);
            end
            Distribution = fnnls(Q,KtS,InitialGuess,TolFun);
            %In some cases, fnnls may return negatives if tolerance is to high
            if any(Distribution < 0)
                %... in those cases continue from current solution
                Distribution = fnnls(Q,KtS,Distribution,1e-20);
            end
    end
    %Store current convergence curve point
    ConvergenceCurve(Iteration) = std(Kernel*Distribution - Signal);
    
    %If hook to axes is given, then plot the current Distribution
    if ~isempty(AxisHandle)
        set(AxisHandle,'YData',Distribution)
        drawnow
    end
    %Update subgradient at current solution
    Subgradient = Subgradient + Kernel'*(Kernel*Distribution - Signal);
    
    
    
    %--------------------------------------------------------------------------
    %Iteration Control
    %--------------------------------------------------------------------------
    if Iteration == 1
        %If at first iteration, thae residual deviation is already below the noise deviation then impose oversmoothing and remain at first iteration
        if NoiseLevelAim  > std(Kernel*Distribution - Signal)
            RegParam = RegParam*2^Counter;
            Counter = Counter + 1;
        else
            %Once the residual deviation is above the treshold, then proceed further with the Bregman iterations
            Iteration  = Iteration +1;
        end
    else
        %For the rest of the Bregman iterations control the condition and stop when fulfilled
        if NoiseLevelAim  > std(Kernel*Distribution - Signal)
            break;
        else
            Iteration  = Iteration +1;
        end
        %If residual deviation starts to diverge, stop
        if DivergenceStop && std(Kernel*CheckDistribution - Signal) < std(Kernel*Distribution - Signal)
            Distribution = CheckDistribution;
            break;
        end
    end
    
end

%Normalize distribution integral
Distribution = Distribution/sum(Distribution)/mean(diff(DistanceAxis));

end

function [Functional,Gradient] = OBIRFunctional(Distribution,RegFunctional,Subgradient)
[FunctionalPart,GradientPart] =  RegFunctional(Distribution);
Functional = FunctionalPart + dot(Distribution,Subgradient);
Gradient = GradientPart + Subgradient;
end