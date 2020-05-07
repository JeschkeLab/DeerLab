function [pass,maxerr] = test(opt)

% Check that fitsignal() works with a form factor (no background)

rng(1)

t = linspace(0,5,100);
r = linspace(2,6,200);
P = dd_gauss(r,[4.5 0.6]);
lam = 0.3;
V = dipolarsignal(t,r,P,lam,'noiselevel',0.01);

[Vfit,Pfit] = fitsignal(V,t,r,'P','none');

% Pass 1: signal is well fitted
pass(1) = all(abs(Vfit - V) < 3e-2);
% Pass 2: distribution is well fitted
pass(2) = all(abs(Pfit - P) < 2e-1);

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
   plot(t,V,'k.',t,Vfit,'r')
   legend('data','fit')
   xlabel('t [\mus]')
   ylabel('V(t)')
   grid on, axis tight, box on
end

end