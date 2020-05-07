function [pass,maxerr] = test(opt)

% Check that fitsignal() fits multiple form factors

rng(1)

t1 = linspace(-0.5,3,500);
t2 = linspace(-0.3,7,300);

r = linspace(2,6,200);
P = dd_gauss(r,[4.5 0.6]);

lam1 = 0.3;
lam2 = 0.4;

V1 = dipolarsignal(t1,r,P,lam1,'noiselevel',0.01);
V2 = dipolarsignal(t2,r,P,lam2,'noiselevel',0.01);

[Vfit,Pfit] = fitsignal({V1,V2},{t1,t2},r,'P','none',@ex_4pdeer);

% Pass 1: signal is well fitted
pass(1) = all(abs(Vfit{1} - V1) < 3e-1);
pass(1) = all(abs(Vfit{2} - V2) < 3e-1);
% Pass 2: distribution is well fitted
pass(2) = all(abs(Pfit - P) < 4e-1);

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
   hold on
   plot(t1,V1,'k.',t1,Vfit{1},'r')
   plot(t2,V2+0.5,'k.',t2,Vfit{2}+0.5,'r')
   legend('data','fit')
   xlabel('t [\mus]')
   ylabel('V(t)')
   grid on, axis tight, box on
end

end