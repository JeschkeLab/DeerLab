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
    GradObj = false;
else
    GradObj = true;
end
validateattributes(RegMatrix,{'numeric'},{'nonempty','2d'},mfilename,'RegMatrix')
validateattributes(RegParam,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'RegParam')
validateattributes(Signal,{'numeric'},{'nonempty'},mfilename,'Signal')
checklengths(Signal,Kernel);

%--------------------------------------------------------------------------
% Parse & Validate Optional Input
%--------------------------------------------------------------------------
%Check if user requested some options via name-value input
[nonNegLSQsolTol,Solver,NonNegConstrained,MaxFunEvals,MaxIter,HuberParam] = parseoptional({'nonNegLSQsolTol','Solver','NonNegConstrained','MaxFunEvals','MaxIter','HuberParam'},varargin);

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
if ~iscolumn(Signal)
    Signal = Signal.';
end

%--------------------------------------------------------------------------
%Memoization
%--------------------------------------------------------------------------

persistent cachedData
if isempty(cachedData)
    cachedData =  java.util.Hashtable;
end

hashKey = datahash({Signal,Kernel,RegMatrix,RegType,RegParam,varargin});
if cachedData.containsKey(hashKey)
    Output = cachedData.get(hashKey);
    [Distribution] = java2mat(Output);
    return
end

%--------------------------------------------------------------------------
%Regularization processing
%--------------------------------------------------------------------------

Dimension = length(Signal);
InitialGuess = zeros(Dimension,1);

%If unconstrained regularization is requested then solve analytically
if ~NonNegConstrained && ~strcmp(Solver,'fmincon')
    Solver = 'analyticalSolution';
end

%If using LSQ-based solvers then precompute the KtK and KtS input arguments
if ~strcmp(Solver,'fmincon')
    KtS = Kernel.'*Signal;
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
end

switch Solver
    
    case 'analyticalSolution'
        PseudoInverse = Q\Kernel.';
        Distribution = PseudoInverse*Signal;
        
    case 'lsqnonneg'
        solverOpts = optimset('Display','off','TolX',nonNegLSQsolTol);
        Distribution = lsqnonneg(Q,KtS,solverOpts);
        
    case 'lsqnonlin'
        RegFunctional = @(x) [Kernel*x - Signal  [RegParam*RegMatrix*x;zeros(size(RegMatrix,2) - size(RegMatrix,1),1)]];
        NonNegConst = zeros(Dimension,1);
        solverOpts=optimoptions(@lsqnonlin,'Display','off',...
            'MaxIter',MaxIter,'MaxFunEvals',MaxFunEvals,...
            'TolFun',nonNegLSQsolTol,'DiffMinChange',1e-8,'DiffMaxChange',0.1);
        [Distribution,~,~,exitflag] = lsqnonlin(RegFunctional,InitialGuess,NonNegConst,[],solverOpts);
        if exitflag == 0
            %... if maxIter exceeded (flag =0) then doube iterations and continue from where it stopped
            solverOpts = optimoptions(solverOpts,'MaxIter',2*MaxIter,'MaxFunEvals',2*MaxFunEvals);
            Distribution = lsqnonlin(RegFunctional,Distribution,NonNegConst,[],solverOpts);
        end
        
    case 'fnnls'
        Distribution = fnnls(Q,KtS,InitialGuess,nonNegLSQsolTol);
        %In some cases, fnnls may return negatives if tolerance is to high
        if any(Distribution < 0)
            %... in those cases continue from current solution
            Distribution = fnnls(Q,KtS,Distribution,1e-20);
        end
    case 'bppnnls'
        Distribution = nnls_bpp(Q,KtS,Q\KtS);
        
    case 'fmincon'
        %Constrained Tikhonov/Total variation/Huber regularization
        if NonNegConstrained
            NonNegConst = zeros(Dimension,1);
        else
            NonNegConst = [];
        end
        if ~strcmp(RegType,'custom')
            RegFunctional = regfunctional(RegType,Signal,RegMatrix,Kernel,RegParam,HuberParam);
        end
        fminconOptions = optimoptions(@fmincon,'SpecifyObjectiveGradient',GradObj,'MaxFunEvals',MaxFunEvals,'Display','off','MaxIter',MaxIter);
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

%Store output result in the cache
cachedData = addcache(cachedData,hashKey,Distribution);

%--------------------------------------------------------------------------
end
