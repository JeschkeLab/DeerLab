%
% OBIR Osher's Bregman-iterated regularization method
%
%   P = OBIR(V,K,r,'type',alpha)
%
%   OBIR of the N-point signal (V) to a M-point distance
%   distribution (P) given a M-point distance axis (r) and NxM point kernel
%   (K). The regularization parameter (alpha) controls the regularization 
%   properties.
%
%   The type of regularization employed in OBIR is set by the 'type'
%   input argument. The regularization models implemented in OBIR are:
%       'tikhonov' -   Tikhonov regularization
%       'tv'       -   Total variation regularization
%       'huber'    -   pseudo-Huber regularization
%
%   P = OBIR(...,'Property',Value)
%   Additional (optional) arguments can be passed as name-value pairs.
% 
%  The name-value pairs to be passed as options can be set in any order.
%
%       'RegOrder' - Order of the regularization operator L (default = 2).
%       'NoiseLevelAim' - Level (standard deviation) of noise at which 
%                         Bregman iterations are to stop.
%       'DivergenceStop'- True/false forces Bregman iterations to stop if
%                         the evolution of the fit's standard deviation 
%                         starts to diverge. 
%       'MaxOuterIter' - Maximal number of Bregman iterations.
%       'Solver' - Minimization solver (default = 'fnnls')
%                      'fmincon' - Non linear constrained minimization
%                      'fnnls' - Fast non-negative least-squares
%       'HuberParam' - Huber parameter used in the 'huber' model (default = 1.35).
%       'TolFun' - Optimizer function tolerance
%       'MaxIter' - Maximum number of optimizer iterations
%       'MaxFunEvals' - Maximum number of optimizer function evaluations
%       'AxisHandle' - Axis handle to plot the state of the distance
%                      distribution at each iteration
%       'Verbose' - Display options for the solvers:
%                    'off' - no information displayed
%                    'final' - display solver exit message
%                    'iter-detailed' - display state of solver at each iteration
%                     See MATLAB doc optimoptions for detailed explanation
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function [P,ConvergenceCurve] = obir(V,K,r,RegType,alpha,varargin)

if nargin<3
    error('At least three inputs (S,K,r) are required.')
end

validateattributes(V,{'numeric'},{'nonempty','vector'},mfilename,'S')
validateattributes(K,{'numeric'},{'nonempty'},mfilename,'K')
V = V(:);
r = r(:);
if length(V)~=size(K,1)
    error('K and signal arguments must fulfill size(K,1)==length(S).')
end
if ~isreal(V)
    error('Input signal cannot be complex.')
end

% Get optional parameters
optionalProperties = {'Verbose','NoiseLevelAim','Solver','MaxIter','TolFun',...
  'MaxFunEvals','DivergenceStop','MaxOuterIter','HuberParam','AxisHandle','RegOrder'};
[Verbose,NoiseLevelAim,Solver,MaxIter,TolFun,MaxFunEvals,DivergenceStop,...
  MaxOuterIter,HuberParam,AxisHandle,RegOrder] ...
    = parseoptional(optionalProperties,varargin);

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
    NoiseLevelAim = noiselevel(V);
else
    validateattributes(NoiseLevelAim,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'NoiseLevelAim')
end

if isa(alpha,'char')
    alpha = selregparam(V,K,RegType,alpha);
else
    validateattributes(alpha,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'RegParam')
end

if nargin<3 || isempty(RegType)
    RegType = 'tikhonov';
else
    validateattributes(RegType,{'char'},{'nonempty'},mfilename,'RegType')
    allowedInput = {'tikhonov','tv','huber'};
    RegType = validatestring(RegType,allowedInput);
end

L = regoperator(length(r),RegOrder);

% Parse & validate optional input
%-------------------------------------------------------------------------------
if isempty(Verbose)
    Verbose = 'off';
else
    validateattributes(Verbose,{'char'},{'nonempty'},mfilename,'Verbose')
end

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

% Preparation
%-------------------------------------------------------------------------------

% Initialize
SizeP = length(r);
Subgradient = zeros(SizeP,1);
Counter = 1;
Iteration = 1;
P = zeros(SizeP,1);

% Osher's Bregman Iterations Algorithm
%-------------------------------------------------------------------------------

Dimension = length(r);
InitialGuess = zeros(Dimension,1);
NonNegConst = zeros(Dimension,1);

while Iteration <= MaxOuterIter
    
    % Store previous iteration distribution
    Pprev = P;
    
    switch Solver
        case 'fmincon'
            % Define current minimization problem
            RegFunctional = regfunctional(RegType,V,L,K,alpha,HuberParam);
            fminconFunctional = @(P)obirfunctional(P,RegFunctional,Subgradient);
            fminconOptions = optimset('GradObj','on','MaxFunEvals',MaxFunEvals,...
              'Display',Verbose,'MaxIter',MaxIter);
            % Run minimzation
            P =  fmincon(fminconFunctional,InitialGuess,[],[],[],[],NonNegConst,[],[],fminconOptions);
        case 'fnnls'
            [KtKreg,KtS] = lsqcomponents(V,K,L,alpha,RegType,HuberParam);
            KtS = KtS - Subgradient;
            P = fnnls(KtKreg,KtS,InitialGuess,TolFun,Verbose);
    end
    % Store current convergence curve point
    ConvergenceCurve(Iteration) = std(K*P - V);
    
    % If hook to axes is given, then plot the current P
    if ~isempty(AxisHandle)
        plot(AxisHandle,ConvergenceCurve)
        drawnow
    end
    % Update subgradient at current solution
    Subgradient = Subgradient + K.'*(K*P - V);
        
    % Iteration control
    %--------------------------------------------------------------------------
    if Iteration == 1
        % If at first iteration, the residual deviation is already below the
        % noise deviation then impose oversmoothing and remain at first iteration
        if NoiseLevelAim  > std(K*P - V)
            alpha = alpha*2^Counter;
            Counter = Counter + 1;
        else
            % Once the residual deviation is above the treshold, then proceed
            % further with the Bregman iterations
            Iteration  = Iteration + 1;
        end
    else
        % For the rest of the Bregman iterations control the condition and stop
        % when fulfilled
        if NoiseLevelAim  > std(K*P - V)
            break;
        else
            Iteration  = Iteration + 1;
        end
        % If residual deviation starts to diverge, stop
        if DivergenceStop && std(K*Pprev - V) < std(K*P - V)
            P = Pprev;
            break;
        end
    end
    
end

% Normalize distribution
if ~all(P==0)
P = P/trapz(r,P);
end

end

function [Functional,Gradient] = obirfunctional(P,RegFunctional,Subgradient)
[FunctionalPart,GradientPart] =  RegFunctional(P);
Functional = FunctionalPart + dot(P,Subgradient);
Gradient = GradientPart + Subgradient;
end
