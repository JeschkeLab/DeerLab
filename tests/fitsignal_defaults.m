function [pass,maxerr] = test(opt)

% Check that fitsignal() works with defaults

rng(1)

t = linspace(0,5,100);
r = linspace(2,6,200);
P = dd_gauss(r,[4.5 0.6]);
B = bg_exp(t,0.2);
lam = 0.3;
V = dipolarsignal(t,r,P,lam,B,'noiselevel',0.01);

[Vfit,Pfit,Bfit,parfit] = fitsignal(V,t,r);
lamfit = parfit.ex;
% Pass 1: signal is well fitted
pass(1) = all(abs(Vfit - V) < 3e-2);
% Pass 2: distribution is well fitted
pass(2) = all(abs(Pfit - P) < 2e-1);
% Pass 3: background is well fitted
pass(3) = all(abs(Bfit - B) < 1e-2);

pass = all(pass);

maxerr = max(abs(Pfit - P));

if opt.Display
   subplot(211)
   plot(r,P,'k',r,Pfit,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
   subplot(212)
   plot(t,V,'k.',t,(1-lam)*B,'k--',t,Vfit,'r',t,(1-lamfit)*Bfit,'r--')
   legend('data','fit')
   xlabel('t [\mus]')
   ylabel('V(t)')
   grid on, axis tight, box on
end

end