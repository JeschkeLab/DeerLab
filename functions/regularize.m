function [Distribution] = regularize(Signal,Kernel,Opts)

if ~iscolumn(Signal)
   Signal = Signal';
end
RegMatrix = getRegMatrix(length(Signal),Opts.RegMatrixOrder);

if isempty(Opts.RegParam)
   RegParam = selectRegParam(Signal,Kernel,RegMatrix);
else
   RegParam = Opts.RegParam; 
end

switch Opts.RegPenalty
    case 'tikhonov'
        RegPenaltyGrad = RegParam^2*(RegMatrix'*RegMatrix);
    case 'tv'
    
end

switch Opts.nonNegLSQsolver
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
end

                     
end
