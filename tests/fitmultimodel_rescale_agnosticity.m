function [pass,maxerr] = test(opt)

% Check that fitmultimodel rescaling does not change the results

rng(1)

t = linspace(0,5,100);
r = linspace(2,6,100);
P = dd_gauss(r,[4 0.6]);
B = bg_exp(t,0.1);
K = dipolarkernel(t,r,0.3,B);

scale = 1e3;
V = K*P + whitegaussnoise(t,0.005);

Pfit1 = fitmultimodel(V*scale,K,r,@dd_gauss,3,'aic');
Pfit2 = fitmultimodel(V,K,r,@dd_gauss,3,'aic','rescale',false);

%Pass 1: distance distribution is well fitted
pass(1) = all(abs(Pfit1 - P) < 2e-1);
pass(2) = all(abs(Pfit2 - Pfit1) < 1e-2);
pass = all(pass);

maxerr = max(abs(Pfit1 - P));

if opt.Display
   plot(r,P,'k',r,Pfit1,'r',r,Pfit2,'b')
   legend('truth','rescaled','normalized')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end