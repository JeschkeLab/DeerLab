function [pass,maxerr] = test(opt)

% Check indifference of fitmultigauss() towards input dimensionality

t = linspace(0,5,100);
r = linspace(1,6,100);
P = dd_gauss(r,[4 0.3]);
K = dipolarkernel(t,r);
S = K*P;

[Pfit1,parfit1] = fitmultimodel(S,K,r,@dd_gauss,2,'aic');
[Pfit2,parfit2] = fitmultimodel(S.',K,r.',@dd_gauss,2,'aic');
[Pfit3,parfit3] = fitmultimodel(S.',K,r,@dd_gauss,2,'aic');
[Pfit4,parfit4] = fitmultimodel(S,K,r.',@dd_gauss,2,'aic');
[Pfit5,parfit5] = fitmultimodel(S,t,r.',@dd_gauss,2,'aic');
[Pfit6,parfit6] = fitmultimodel(S,t.',r,@dd_gauss,2,'aic');
[Pfit7,parfit7] = fitmultimodel(S,t.',r.',@dd_gauss,2,'aic');

% Pass 1: all distributions are equal
pass(1) = isequal(Pfit1,Pfit2,Pfit3,Pfit4,Pfit5,Pfit6,Pfit7);
% Pass 2: all distributions are columns vectors
pass(2) = iscolumn(Pfit1) & iscolumn(Pfit2) & iscolumn(Pfit3) & iscolumn(Pfit4) & iscolumn(Pfit5) & iscolumn(Pfit6)& iscolumn(Pfit7);
% Pass 3: all fit parameters are row vectors
pass(3) = ~iscolumn(parfit1) & ~iscolumn(parfit2) & ~iscolumn(parfit3) & ~iscolumn(parfit4) & ~iscolumn(parfit5) & ~iscolumn(parfit6) & ~iscolumn(parfit7);

pass = all(pass);

maxerr = NaN;
 

if opt.Display
   plot(r,P,'k',r,Pfit1,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end