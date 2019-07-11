function Distribution = OBIR(Signal,Kernel,Method,NoiseLevelAim,Opts)
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
%
% L. Fabregas, 2018

%--------------------------------------------------------------------------
% Preparation
%--------------------------------------------------------------------------

%Initialize
SizeDistribution = length(Signal);
Subgradient = zeros(SizeDistribution,1);
Counter = 1;
Iteration = 1;
Signal = Signal';

Distribution = zeros(SizeDistribution,1);

%--------------------------------------------------------------------------
% Osher's Bregman Iterations Algorithm
%--------------------------------------------------------------------------

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

while Iteration <= 200
    CheckDistribution = Distribution;    
    
    InitialGuess = zeros(Dimension,1);
    NonNegConst = zeros(Dimension,1);
    RegFunctional = getRegFunctional(Method,Signal,RegMatrix,Kernel,RegParam);
    BregmanDist = @(Distribution)dot(Distribution,Subgradient);
    OBIRFunctional = @(Distribution) (RegFunctional(Distribution) + BregmanDist(Distribution));
    fminconOptions = optimset('GradObj','off','MaxFunEvals',200000,'Display','off','MaxIter',200000);
    Distribution =  fmincon(OBIRFunctional,InitialGuess,[],[],[],[],NonNegConst,[],[],fminconOptions);
    
    %Update subgradient at current solution
    Subgradient = Subgradient + Kernel'*(Kernel*Distribution - Signal);
    
    Distribution = Distribution/sum(Distribution);
        ResidualDeviation(Iteration) = std(Kernel*Distribution - Signal);
        figure(4124),clf,subplot(1,2,1),plot(Distribution/sum(Distribution),'k')
        subplot(122),plot(1:Iteration,ResidualDeviation,'ko--');hold on,plot(1:Iteration+1,ones(Iteration+1,1)*NoiseLevelAim,'r--');
        xlabel('Iteration'),ylabel('Residual Std. Dev.'),ylims = ylim; ylim([0.95*ylims(1) ylims(2)]),drawnow,
    
    %--------------------------------------------------------------------------
    %Iteration Control
    %--------------------------------------------------------------------------
    if Iteration == 1
        %If at first iteration, the residual deviation is already below the noise deviation then impose oversmoothing and remain at first iteration
        if NoiseLevelAim  > std(Kernel*Distribution - Signal)
            RegParam = RegParam*2^Counter;
            Counter = Counter + 1;
        else
            %Once the residual deviation is above the treshold, then proceed further with the Bregman iterations
            Iteration  = Iteration +1;
        end
    else
        %For the rest of the Bregman iterations then control the condition and stop when fulfilled
        if NoiseLevelAim  > std(Kernel*Distribution - Signal)
            break;
        else
            Iteration  = Iteration +1;
        end
        %If residual deviation starts to diverge, stop
        if std(Kernel*CheckDistribution - Signal) < std(Kernel*Distribution - Signal)
            Distribution = CheckDistribution;
            break;
        end
    end
    
end

end
