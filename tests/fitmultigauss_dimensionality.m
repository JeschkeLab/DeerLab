function [pass,maxerr] = test(opt)

% Check indifference of fitmultigauss() towards input dimensionality

t = linspace(0,5,100);
r = linspace(1,6,100);
P = rd_onegaussian(r,[4 0.3]);
K = dipolarkernel(t,r);
S = K*P;

[Pfit1,parfit1] = fitmultigauss(S,K,r,3,'aic');
[Pfit2,parfit2] = fitmultigauss(S.',K,r.',3,'aic');
[Pfit3,parfit3] = fitmultigauss(S.',K,r,3,'aic');
[Pfit4,parfit4] = fitmultigauss(S,K,r.',3,'aic');

% Pass 1: all distributions are equal
pass(1) = isequal(Pfit1,Pfit2,Pfit3,Pfit4);
% Pass 1: all distributions are columns vectors
pass(2) = iscolumn(Pfit1) & iscolumn(Pfit2) & iscolumn(Pfit3) & iscolumn(Pfit4);
% Pass 1: all fit parameters are row vectors
pass(3) = ~iscolumn(parfit1) & ~iscolumn(parfit2) & ~iscolumn(parfit3) & ~iscolumn(parfit4);

pass = all(pass);

maxerr = max(abs(Pfit1 - Pfit2));
 

if opt.Display
   plot(r,P,'k',r,Pfit1,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end