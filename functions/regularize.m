function Distribution = regularize(Signal,Kernel,RegMatrix,RegType,RegParam,varargin)

%--------------------------------------------------------------------------
% Parse & Validate Required Input
%--------------------------------------------------------------------------
if nargin<5
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
    validatestring(RegType,allowedInput);
end

if ~isreal(Signal)
    Signal = real(Signal);
end

if strcmp(RegType,'custom')
    GradObj = 'off';
else
    GradObj = 'on';
end
validateattributes(RegMatrix,{'numeric'},{'nonempty','2d'},mfilename,'RegMatrix')
validateattributes(RegParam,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'RegParam')
validateattributes(Signal,{'numeric'},{'nonempty'},mfilename,'Signal')
checklengths(Signal,Kernel);

%--------------------------------------------------------------------------
% Parse & Validate Optional Input
%--------------------------------------------------------------------------
%Check if user requested some options via name-value input
[nonNegLSQsolTol,Solver,NonNegConstrained,MaxFunEvals,MaxIter] = parseoptional({'nonNegLSQsolTol','Solver','NonNegConstrained','MaxFunEvals','MaxIter'},varargin);

if isempty(nonNegLSQsolTol)
    nonNegLSQsolTol = 1e-9;
else
    validateattributes(nonNegLSQsolTol,{'numeric'},{'scalar','nonempty','nonnegative'},'regularize','nonNegLSQsolTol')
end
if isempty(Solver)
    Solver = 'fmincon';
else
    validateattributes(Solver,{'char'},{'nonempty'})
    allowedInput = {'fnnls','lsqnonneg','bppnnls','fmincon','cvx','lsqnonlin'};
    validatestring(Solver,allowedInput);
end

if isempty(MaxIter)
    MaxIter = 20000000;
else
    validateattributes(MaxIter,{'numeric'},{'scalar','nonempty'},mfilename,'MaxIter')
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
if ~iscolumn(Signal)
    Signal = Signal';
end

checkSolverCompatibility(Solver,RegType);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Regularization processing
%--------------------------------------------------------------------------

Dimension = length(Signal);
InitialGuess = zeros(Dimension,1);

%If unconstrained Tikh regularization is requested then solve analytically
if ~NonNegConstrained && strcmp(RegType,'tikhonov')
    Solver = 'analyticaluTikhonov';
end

switch Solver
    
    case 'analyticaluTikhonov'
        %Analytical solution of unconstrained Tikhonov problem
        PseudoInverse = ((Kernel.'*Kernel) + RegParam^2*(RegMatrix.'*RegMatrix))\(Kernel.');
        Distribution = PseudoInverse*Signal;
        
    case 'lsqnonneg'
        %Constrained Tikhonov regularization
        solverOpts = optimset('Display','off','TolX',nonNegLSQsolTol);
        Q = (Kernel'*Kernel) + RegParam^2*(RegMatrix'*RegMatrix);
        Distribution = lsqnonneg(Q,Kernel'*Signal,solverOpts);
    case 'lsqnonlin'
        RegFunctional = @(x) [Kernel*x - Signal  [RegParam*RegMatrix*x;0;0]];
        NonNegConst = zeros(Dimension,1);
        solverOpts=optimoptions(@lsqnonlin,'Display','off',...
            'MaxIter',MaxIter,'MaxFunEvals',MaxFunEvals,...
            'TolFun',nonNegLSQsolTol,'DiffMinChange',1e-8,'DiffMaxChange',0.1);
        [Distribution,~,~,exitflag] = lsqnonlin(RegFunctional,InitialGuess,NonNegConst,[],solverOpts);
        if exitflag == 0
            %... if maxIter exceeded (flag =0) then doube iterations and continue from where it stopped
            solverOpts=optimoptions(solverOpts,'MaxIter',2*MaxIter,'MaxFunEvals',2*MaxFunEvals);
            Distribution = lsqnonlin(RegFunctional,Distribution,NonNegConst,[],solverOpts);
        end
    case 'fnnls'
        %Constrained Tikhonov regularization
        Q = (Kernel'*Kernel) + RegParam^2*(RegMatrix'*RegMatrix);
        Distribution = fnnls(Q,Kernel'*Signal,InitialGuess,nonNegLSQsolTol);
        
    case 'bppnnls'
        %Constrained Tikhonov regularization
        Q = (Kernel'*Kernel) + RegParam^2*(RegMatrix'*RegMatrix);
        KtS = Kernel'*Signal;
        Distribution = nnls_bpp(Q,KtS,Q\KtS);
        
    case 'fmincon'
        %Constrained Tikhonov/Total variation/Huber regularization
        if NonNegConstrained
            NonNegConst = zeros(Dimension,1);
        else
            NonNegConst = [];
        end
        if ~strcmp(RegType,'custom')
            RegFunctional = regfunctional(RegType,Signal,RegMatrix,Kernel,RegParam);
        end
        fminconOptions = optimoptions(@fmincon,'SpecifyObjectiveGradient',true,'MaxFunEvals',MaxFunEvals,'Display','off','MaxIter',MaxIter);
        [Distribution,~,exitflag] =  fmincon(RegFunctional,InitialGuess,[],[],[],[],NonNegConst,[],@unityconstraint,fminconOptions);
        %Check how optimization exited...
        if exitflag == 0
            %... if maxIter exceeded (flag =0) then doube iterations and continue from where it stopped
            fminconOptions = optimoptions(fminconOptions,'MaxIter',2*MaxIter,'MaxFunEvals',2*MaxFunEvals);
            Distribution  = fmincon(ModelCost,Distribution,[],[],[],[],NonNegConst,[],@unityconstraint,fminconOptions);
        end
    case 'cvx'
        %Constrained Tikhonov/Total variation/Huber regularization
        Distribution = cvxregsolver(RegType,Signal,RegMatrix,Kernel,RegParam,NonNegConstrained);
end

%Normalize distribution integral
Distribution = Distribution/sum(Distribution);
%--------------------------------------------------------------------------
end
