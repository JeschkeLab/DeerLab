% 
% FITREGMODEL Fits a distance distribution to one (or several) signals
%            by optimization of a regularization functional model.
%
%   P = FITREGMODEL(S,K,r,L,regtype,alpha)
%   Regularization of the N-point signal (S) to a M-point distance
%   distribution (P) given a M-point distance axis (r) and NxM point kernel
%   (K). The (M-2)xM point regularization matrix (L) and regularization
%   parameter (alpha) control the regularization properties.
%
%   The type of regularization employed in FITREGMODEL is set by the regtype
%   input argument. The regularization models implemented in FITREGMODEL are:
%          'tikhonov' -   Tikhonov regularization
%          'tv'       -   Total variation regularization
%          'huber'    -   pseudo-Huber regularization
%
%   P = FITREGMODEL({S1,S2,...},{K1,K2,...},r,L,regtype,alpha)
%   Passing multiple signals/kernels enables global fitting of the
%   regularization model to a single distribution. The global fit weights
%   are automatically computed according to their contribution to ill-posedness.
%
%   P = FITREGMODEL(...,'Property',Values)
%   Additional (optional) arguments can be passed as property-value pairs.
%
% The properties to be passed as options can be set in any order. 
%
%   'Solver' - Solver to be used to solve the minimization problems
%                      'fnnls' - Fast non-negative least-squares
%                      'lsqnonneg' - Non-negative least-squares 
%                      'fmincon' - Non-linear constrained minimization
%                      'bppnnls' -  Block principal pivoting non-negative least-squares solver
%
%   'NonNegConstrained' - Enable/disable non-negativity constraint (true/false)
%
%   'HuberParam' - Huber parameter used in the 'huber' model (default = 1.35).
%
%   'GlobalWeights' - Array of weighting coefficients for the individual signals in
%                     global fitting regularization.
%
%   'TolFun' - Optimizer function tolerance
%
%   'MaxIter' - Maximum number of optimizer iterations
%
%   'MaxFunEvals' - Maximum number of optimizer function evaluations   
%
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function P = fitregmodel(S,K,r,RegMatrix,RegType,RegParam,varargin)


%Turn off warnings to avoid ill-conditioned warnings 
warning('off','MATLAB:nearlySingularMatrix')

%--------------------------------------------------------------------------
% Parse & Validate Required Input
%--------------------------------------------------------------------------
if nargin<6
    error('Not enough input arguments.')
end
if nargin<4 || isempty(RegType)
    RegType = 'tikhonov';
elseif isa(RegType,'function_handle')
    RegFunctional = RegType;
    RegType = 'custom';
else
    validateattributes(RegType,{'char'},{'nonempty'})
    allowedInput = {'tikhonov','tv','huber'};
    RegType = validatestring(RegType,allowedInput);
end

if strcmp(RegType,'custom')
    GradObj = false;
else
    GradObj = true;
end
validateattributes(RegMatrix,{'numeric'},{'nonempty','2d'},mfilename,'RegMatrix')
validateattributes(RegParam,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'RegParam')
validateattributes(r,{'numeric'},{'nonempty','increasing','nonnegative'},mfilename,'r')
if numel(unique(round(diff(r),6)))~=1
    error('Distance axis must be a monotonically increasing vector.')
end

%--------------------------------------------------------------------------
% Parse & Validate Optional Input
%--------------------------------------------------------------------------
%Check if user requested some options via name-value input
[TolFun,Solver,NonNegConstrained,MaxFunEvals,MaxIter,HuberParam,GlobalWeights] = parseoptional({'TolFun','Solver','NonNegConstrained','MaxFunEvals','MaxIter','HuberParam','GlobalWeights'},varargin);

if isempty(TolFun)
    TolFun = 1e-9;
else
    validateattributes(TolFun,{'numeric'},{'scalar','nonempty','nonnegative'},'regularize','nonNegLSQsolTol')
end
if isempty(Solver)
    Solver = 'fnnls';
else
    validateattributes(Solver,{'char'},{'nonempty'})
    allowedInput = {'analytical','fnnls','lsqnonneg','bppnnls','fmincon'};
    Solver = validatestring(Solver,allowedInput);
end

if isempty(MaxIter)
    MaxIter = 20000000;
else
    validateattributes(MaxIter,{'numeric'},{'scalar','nonempty'},mfilename,'MaxIter')
end

if isempty(HuberParam)
    HuberParam = 1.35;
else
    validateattributes(HuberParam,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'MaxFunEvals')
end

if isempty(MaxFunEvals)
    MaxFunEvals = 2000000;
else
    validateattributes(MaxFunEvals,{'numeric'},{'scalar','nonempty'},mfilename,'MaxFunEvals')
end

if isempty(NonNegConstrained)
    NonNegConstrained = true;
else
    validateattributes(NonNegConstrained,{'logical'},{'nonempty'},'regularize','NonNegConstrained')
end
if ~iscell(S)
    S = {S};
end
if ~iscell(K)
    K = {K};
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
if length(K)~=length(S)
    error('The number of kernels and signals must be equal.')
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
end

%--------------------------------------------------------------------------
%Regularization processing
%--------------------------------------------------------------------------

Dimension = length(RegMatrix);
InitialGuess = zeros(Dimension,1);

%Convert distance axis to nanoseconds if givne in Angstrom
if ~isnanometer(r)
   r = r/10; 
end
dr = mean(diff(r));

%If unconstrained regularization is requested then solve analytically
if ~NonNegConstrained && ~strcmp(Solver,'fmincon')
    Solver = 'analytical';
end

%If using LSQ-based solvers then precompute the KtK and KtS input arguments
if ~strcmp(Solver,'fmincon')
    [Q,KtS,weights] =  lsqcomponents(S,K,RegMatrix,RegParam,RegType,HuberParam,GlobalWeights);
end

%Solve the regularization functional minimization problem
switch lower(Solver)
    
    case 'analytical'
        P = zeros(Dimension,1);
        for i=1:length(S)
        PseudoInverse = Q\K{i}.';
        P = P + weights(i)*PseudoInverse*S{i};
        end
    case 'lsqnonneg'
        solverOpts = optimset('Display','off','TolX',TolFun);
        P = lsqnonneg(Q,KtS,solverOpts);

    case 'fnnls'
        P = fnnls(Q,KtS,InitialGuess,TolFun);
        %In some cases, fnnls may return negatives if tolerance is to high
        if any(P < 0)
            %... in those cases continue from current solution
            P = fnnls(Q,KtS,P,1e-20);
        end
    case 'bppnnls'
        P = nnls_bpp(Q,KtS,Q\KtS);
        
    case 'fmincon'
        %Constrained Tikhonov/Total variation/Huber regularization
        if NonNegConstrained
            NonNegConst = zeros(Dimension,1);
        else
            NonNegConst = [];
        end
        if ~strcmp(RegType,'custom')
            RegFunctional = regfunctional(RegType,S,RegMatrix,K,RegParam,HuberParam);
        end
        constraint = @(x)unityconstraint(x,dr);
        fminconOptions = optimoptions(@fmincon,'SpecifyObjectiveGradient',GradObj,'MaxFunEvals',MaxFunEvals,'Display','off','MaxIter',MaxIter);
        [P,~,exitflag] =  fmincon(RegFunctional,InitialGuess,[],[],[],[],NonNegConst,[],constraint,fminconOptions);
        %Check how optimization exited...
        if exitflag == 0
            %... if maxIter exceeded (flag =0) then doube iterations and continue from where it stopped
            fminconOptions = optimoptions(fminconOptions,'MaxIter',2*MaxIter,'MaxFunEvals',2*MaxFunEvals);
            P  = fmincon(ModelCost,P,[],[],[],[],NonNegConst,[],constraint,fminconOptions);
        end
end

%Normalize distribution integral
P = P/sum(P)/dr;

%Turn warnings back on
warning('on','MATLAB:nearlySingularMatrix')

end
