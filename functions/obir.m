%
% OBIR Osher's Bregman-iterated regularization method
%
%   P = OBIR(S,K,r,'type',alpha)
%   OBIR of the N-point signal (S) to a M-point distance
%   distribution (P) given a M-point distance axis (r) and NxM point kernel
%   (K). The regularization parameter (alpha) controls the regularization 
%   properties.
%
%   The type of regularization employed in OBIR is set by the 'type'
%   input argument. The regularization models implemented in OBIR are:
%          'tikhonov' -   Tikhonov regularization
%          'tv'       -   Total variation regularization
%          'huber'    -   pseudo-Huber regularization
%
%   P = OBIR(...,'Property',Value)
%   Additional (optional) arguments can be passed as property-value pairs.
% 
%  The property-value pairs to be passed as options can be set in any order.
%
%       'RegOrder' - Order of the regularization operator L (default = 2).
%
%       'NoiseLevelAim' - Level (standard deviation) of noise at which 
%                         Bregman iterations are to stop.
%
%       'DivergenceStop'- True/false forces Bregman iterations to stop if
%                         the evolution of the fit's standard deviation 
%                         starts to diverge. 
%
%       'MaxOuterIter' - Maximal number of Bregman iterations.
%
%       'Solver' - Minimization solver (default = 'fnnls')
%                      'fmincon' - Non linear constrained minimization
%                      'fnnls' - Fast non-negative least-squares
%
%       'TolFun' - Optimizer function tolerance
%
%       'MaxIter' - Maximum number of optimizer iterations
%
%       'MaxFunEvals' - Maximum number of optimizer function evaluations   
%
%       'AxisHandle' - Axis handle to plot the state of the distance
%                      distribution at each iteration     
%

% This file is a part of DeerAnalysis. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.


function [P,ConvergenceCurve] = obir(S,K,r,RegType,alpha,varargin)

if ~iscolumn(S)
    S = S';
end

%Get optional parameters
[NoiseLevelAim,Solver,MaxIter,TolFun,MaxFunEvals,DivergenceStop,MaxOuterIter,HuberParam,AxisHandle,RegOrder] ...
    = parseoptional({'NoiseLevelAim','Solver','MaxIter','TolFun','MaxFunEvals','DivergenceStop','MaxOuterIter','HuberParam','AxisHandle','RegOrder'},varargin);

if isempty(RegOrder)
    RegOrder = 2;
else
    validateattributes(RegOrder,{'numeric'},{'scalar','nonnegative'})
end

if isempty(HuberParam)
    HuberParam = 1.35;
else
    validateattributes(HuberParam,{'numeric'},{'scalar','nonnegative'},mfilename,'HuberParam')
end

if isempty(Solver)
    Solver = 'fnnls';
else
    validateattributes(Solver,{'char'},{'nonempty'},mfilename,'Solver')
    validInputs = {'fnnls','fmincon'};
    Solver = validatestring(Solver,validInputs);
end

if isempty(NoiseLevelAim)
    NoiseLevelAim = noiselevel(S);
else
    validateattributes(NoiseLevelAim,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'NoiseLevelAim')
end

if isa(alpha,'char')
    alpha = selregparam(S,K,r,RegType,alpha);
else
    validateattributes(alpha,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'RegParam')
end
validateattributes(S,{'numeric'},{'nonempty'},mfilename,'S')
validateattributes(K,{'numeric'},{'nonempty'},mfilename,'K')
checklengths(S,K);
if ~isreal(S)
    error('Input signal cannot be complex.')
end

if nargin<3 || isempty(RegType)
    RegType = 'tikhonov';
else
    validateattributes(RegType,{'char'},{'nonempty'},mfilename,'RegType')
    allowedInput = {'tikhonov','tv','huber'};
    RegType = validatestring(RegType,allowedInput);
end

%Convert distance axis to nanoseconds if givne in Angstrom
if ~isnanometer(r)
   r = r/10; 
end

L = regoperator(length(r),RegOrder);

%--------------------------------------------------------------------------
% Parse & Validate Optional Input
%--------------------------------------------------------------------------

if isempty(TolFun)
    TolFun = 1e-10;
else
    validateattributes(TolFun,{'numeric'},{'scalar','nonnegative'},mfilename,'TolFun')
end

if isempty(MaxOuterIter)
    MaxOuterIter = 1000;
else
    validateattributes(MaxOuterIter,{'numeric'},{'scalar','nonempty'},mfilename,'MaxOuterIter')
end

if isempty(MaxIter)
    MaxIter = 500000;
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
SizeP = length(r);
Subgradient = zeros(SizeP,1);
Counter = 1;
Iteration = 1;
P = zeros(SizeP,1);

%--------------------------------------------------------------------------
% Osher's Bregman Iterations Algorithm
%--------------------------------------------------------------------------

Dimension = length(r);
InitialGuess = zeros(Dimension,1);
NonNegConst = zeros(Dimension,1);

while Iteration <= MaxOuterIter
    
    %Store privous iteration distribution
    CheckP = P;
    
    switch Solver
        case 'fmincon'
            %Define current minimization problem
            RegFunctional = regfunctional(RegType,S,L,K,alpha,HuberParam);
            fminconFunctional = @(P)obirfunctional(P,RegFunctional,Subgradient);
            fminconOptions = optimset('GradObj','on','MaxFunEvals',MaxFunEvals,'Display','off','MaxIter',MaxIter);
            %Run minimzation
            P =  fmincon(fminconFunctional,InitialGuess,[],[],[],[],NonNegConst,[],[],fminconOptions);
        case 'fnnls'
            [Q,KtS] = lsqcomponents(S,r,K,L,alpha,RegType,HuberParam);
            KtS = KtS - Subgradient;
            P = fnnls(Q,KtS,InitialGuess,TolFun);
    end
    %Store current convergence curve point
    ConvergenceCurve(Iteration) = std(K*P - S);
    
    %If hook to axes is given, then plot the current P
    if ~isempty(AxisHandle)
        plot(AxisHandle,ConvergenceCurve),drawnow
    end
    %Update subgradient at current solution
    Subgradient = Subgradient + K'*(K*P - S);
    
    
    
    %--------------------------------------------------------------------------
    %Iteration Control
    %--------------------------------------------------------------------------
    if Iteration == 1
        %If at first iteration, thae residual deviation is already below the noise deviation then impose oversmoothing and remain at first iteration
        if NoiseLevelAim  > std(K*P - S)
            alpha = alpha*2^Counter;
            Counter = Counter + 1;
        else
            %Once the residual deviation is above the treshold, then proceed further with the Bregman iterations
            Iteration  = Iteration +1;
        end
    else
        %For the rest of the Bregman iterations control the condition and stop when fulfilled
        if NoiseLevelAim  > std(K*P - S)
            break;
        else
            Iteration  = Iteration +1;
        end
        %If residual deviation starts to diverge, stop
        if DivergenceStop && std(K*CheckP - S) < std(K*P - S)
            P = CheckP;
            break;
        end
    end
    
end

%Normalize distribution integral
P = P/sum(P)/mean(diff(r));

end

function [Functional,Gradient] = obirfunctional(P,RegFunctional,Subgradient)
[FunctionalPart,GradientPart] =  RegFunctional(P);
Functional = FunctionalPart + dot(P,Subgradient);
Gradient = GradientPart + Subgradient;
end