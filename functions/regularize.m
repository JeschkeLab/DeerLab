function [Distribution] = regularize(Signal,Kernel,Method,Opts)

if nargin<3 || isempty(Method)
    Method = 'tikhonov';
end
if nargin<4 || isempty(Opts)
    Opts = DAoptions;
else
    if ~isa(Opts,'DAoptions')
        error('First argument must a valid DAoptions class object')
    end
end

checkSolverCompatibility(Opts,Method);

if ~iscolumn(Signal)
   Signal = Signal';
end
Dimension = length(Signal);
RegMatrix = getRegMatrix(Dimension,Opts.RegMatrixOrder);

if isempty(Opts.RegParam)
   RegParam = selectRegParam(Signal,Kernel,RegMatrix);
else
   RegParam = Opts.RegParam; 
end

switch Method
    case 'tikhonov'
        RegPenaltyGrad = RegParam^2*(RegMatrix'*RegMatrix);
    case 'tv'
    
end

switch Opts.Solver
  case 'lsqnonneg'
    solverOpts = optimset('Display','off','TolX',Opts.nonNegLSQsolTol);
    Q = (Kernel'*Kernel) + RegPenaltyGrad;
    Distribution = lsqnonneg(Q,Kernel'*Signal,solverOpts);
    
  case 'fnnls'
    Q = (Kernel'*Kernel) + RegPenaltyGrad;
    Distribution = fnnls(Q,Kernel'*Signal,[],Opts.nonNegLSQsolTol);
    
  case 'bppnnls'
    Q = (Kernel'*Kernel) + RegPenaltyGrad;
    KtS = Kernel'*Signal;
    Distribution = nnls_bpp(Q,KtS,Q\KtS);
    
  case 'fmincon'
      InitialGuess = zeros(Dimension,1);
      NonNegConst = zeros(Dimension,1);
      RegFunctional = getRegFunctional(Method,Signal,RegMatrix,Kernel,RegParam);
      fminconOptions = optimset('GradObj','on','MaxFunEvals',200000,'Display','off','MaxIter',200000);
      Distribution =  fmincon(RegFunctional,InitialGuess,[],[],[],[],NonNegConst,[],[],fminconOptions);

end



                     
end
