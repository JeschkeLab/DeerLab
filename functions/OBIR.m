function [Distribution,ConvergenceCurve] = obir(Signal,Kernel,RegType,RegMatrix,RegParam,NoiseLevelAim,varargin)
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
validateattributes(NoiseLevelAim,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'NoiseLevelAim')
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
%Check if user requested some options via name-value input
[MaxOuterIter,MaxFunEvals,MaxIter,DivergenceStop] = parseOptional({'MaxOuterIter','MaxFunEvals','MaxIter','DivergenceStop'},varargin);

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
    
    %Define current minimization problem
    RegFunctional = getRegFunctional(RegType,Signal,RegMatrix,Kernel,RegParam);
    fminconFunctional = @(Distribution)OBIRFunctional(Distribution,RegFunctional,Subgradient);
    fminconOptions = optimset('GradObj','on','MaxFunEvals',MaxFunEvals,'Display','off','MaxIter',MaxIter);
    
    %Run minimzation
    Distribution =  fmincon(fminconFunctional,InitialGuess,[],[],[],[],NonNegConst,[],@unitIntConstraint,fminconOptions);
    
    %Store current convergence curve point
    ConvergenceCurve(Iteration) = std(Kernel*Distribution - Signal);
    
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

end

function [Functional,Gradient] = OBIRFunctional(Distribution,RegFunctional,Subgradient)
[FunctionalPart,GradientPart] =  RegFunctional(Distribution);
Functional = FunctionalPart + dot(Distribution,Subgradient);
Gradient = GradientPart + Subgradient;
end