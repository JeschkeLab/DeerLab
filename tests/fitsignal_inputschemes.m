function [pass,maxerr] = test(opt)

% Check that the different input schemes are equivalent

rng(1)

t = linspace(0,5,100);
r = linspace(2,6,30);
P = dd_gauss(r,[4.5 0.6]);
lam = 0.3;
B = bg_hom3d(t,50,lam);
V = dipolarsignal(t,r,P,lam,B,'noiselevel',0.01);

[~,Pfit1] = fitsignal(V,t,r);
[~,Pfit2] = fitsignal(V,t,r,'P');
[~,Pfit3] = fitsignal(V,t,r,'P',@bg_hom3d);
[~,Pfit4] = fitsignal(V,t,r,'P',@bg_hom3d,@ex_4pdeer);

% Pass 1: al input schemes yield the same results
pass = isequal(Pfit1,Pfit2,Pfit3,Pfit4);

maxerr = NaN;

if opt.Display
   plot(r,P,'k',r,Pfit1,r,Pfit2,r,Pfit3,r,Pfit4)
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end