function Distribution = regularize(Signal,Kernel,RegType,RegParam,varargin)

%--------------------------------------------------------------------------
%Input parsing & validation
%--------------------------------------------------------------------------
if nargin<3 || isempty(RegType)
    RegType = 'tikhonov';
elseif isa(RegType,'function_handle')
    RegFunctional = RegType;
    RegType = 'custom';
else
    validateattributes(RegType,{'char'},{'nonempty'})
    allowedInput = {'tikhonov','tv','huber'};
    validatestring(RegType,allowedInput);
end
if strcmp(RegType,'custom')
    GradObj = 'off';
else
    GradObj = 'on';
end
%Check if user requested some options via name-value input
[RegMatrixOrder,nonNegLSQsolTol,Solver,NonNegConstrained,MaxFunEvals,MaxIter] = parseOptional({'RegMatrixOrder','nonNegLSQsolTol','Solver','NonNegConstrained','MaxFunEvals','MaxIter'},varargin);

if isempty(RegMatrixOrder)
    RegMatrixOrder = 2;
else
    validateattributes(RegMatrixOrder,{'numeric'},{'nonempty','nonnegative'},'regularize','RegMatrixOrder')
end
if isempty(nonNegLSQsolTol)
    nonNegLSQsolTol = 1e-9;
else
    validateattributes(nonNegLSQsolTol,{'numeric'},{'scalar','nonempty','nonnegative'},'regularize','nonNegLSQsolTol')
end
if isempty(Solver)
    Solver = 'fmincon';
else
    validateattributes(Solver,{'char'},{'nonempty'})
    allowedInput = {'fnnls','lsqnonneg','bppnnls','fmincon','cvx'};
    validatestring(Solver,allowedInput);
end
if isempty(NonNegConstrained)
    NonNegConstrained = true;
else
    validateattributes(NonNegConstrained,{'logical'},{'nonempty'},'regularize','NonNegConstrained')
end
validateattributes(RegParam,{'numeric'},{'scalar','nonempty','nonnegative'},'regularize','RegParam')
validateattributes(Signal,{'numeric'},{'nonempty'},'regularize','Signal')
checklengths(Signal,Kernel);
checkSolverCompatibility(Solver,RegType);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Regularization processing
%--------------------------------------------------------------------------

%Get regularization matrix
Dimension = length(Signal);
RegMatrix = getRegMatrix(Dimension,RegMatrixOrder);

%If unconstrained Tikh regularization is requested then solve analytically
if ~NonNegConstrained && strcmp(RegType,'tikhonov')
    Solver = 'analyticaluTikhonov';
end

switch Solver
    
    case 'analyticaluTikhonov'
        %Analytical solution of unconstrained Tikhonov problem
        PseudoInverse = (Kernel'*Kernel) + RegParam^2*(RegMatrix'*RegMatrix)\Kernel';
        Distribution = PseudoInverse*Signal;
        
    case 'lsqnonneg'
        %Constrained Tikhonov regularization
        solverOpts = optimset('Display','off','TolX',nonNegLSQsolTol);
        Q = (Kernel'*Kernel) + RegParam^2*(RegMatrix'*RegMatrix);
        Distribution = lsqnonneg(Q,Kernel'*Signal,solverOpts);
        
    case 'fnnls'
        %Constrained Tikhonov regularization
        Q = (Kernel'*Kernel) + RegParam^2*(RegMatrix'*RegMatrix);
        Distribution = fnnls(Q,Kernel'*Signal,[],nonNegLSQsolTol);
        
    case 'bppnnls'
        %Constrained Tikhonov regularization
        Q = (Kernel'*Kernel) + RegParam^2*(RegMatrix'*RegMatrix);
        KtS = Kernel'*Signal;
        Distribution = nnls_bpp(Q,KtS,Q\KtS);
        
    case 'fmincon'
        %Constrained Tikhonov/Total variation/Huber regularization
        InitialGuess = zeros(Dimension,1);
        if NonNegConstrained
            NonNegConst = zeros(Dimension,1);
        else
            NonNegConst = [];
        end
        if ~strcmp(RegType,'custom')
            RegFunctional = getRegFunctional(RegType,Signal,RegMatrix,Kernel,RegParam);
        end
        fminconOptions = optimset('GradObj',GradObj,'MaxFunEvals',MaxFunEvals,'Display','off','MaxIter',MaxIter);
        Distribution =  fmincon(RegFunctional,InitialGuess,[],[],[],[],NonNegConst,[],[],fminconOptions);
        
    case 'cvx'
        %Constrained Tikhonov/Total variation/Huber regularization
        Distribution = cvxregsolver(RegType,Signal,RegMatrix,Kernel,RegParam,NonNegConstrained);
end

%Normalize distribution integral
Distribution = Distribution/sum(Distribution);
%--------------------------------------------------------------------------
end
