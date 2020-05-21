function [pass,maxerr] = test(opt)

% Check that fitregmodel rescaling does not change the results

rng(1)

t = linspace(0,5,100);
r = linspace(2,6,100);
P = dd_gauss(r,[4 0.6]);
B = bg_exp(t,0.1);
K = dipolarkernel(t,r,0.3,B);

scale = 1e9;
V = K*P + whitegaussnoise(t,0.005);

Pfit1 = fitregmodel(V*scale,K,r,'tikh','aic','normP',false);
Pfit2 = fitregmodel(V,K,r,'tikh','aic');

%Pass 1-3: distance distributions are well fitted
pass(1) = all(abs(Pfit1 - scale*P) < scale*2e-1);
pass(2) = all(abs(Pfit2 - P) < 2e-1);
pass(3) = all(abs(scale*Pfit2 - Pfit1) < scale*1e-3);

pass = all(pass);

maxerr = max(abs(Pfit1 - scale*P));

if opt.Display
   plot(r,P,'k',r,Pfit1,'r',r,Pfit2,'b')
   legend('truth','rescaled','normalized')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end