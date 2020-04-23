function [pass,maxerr] = test(opt)

% Check that fitsignal() works with a dipolar evolution function

rng(1)

t = linspace(0,5,100);
r = linspace(2,6,70);
P = dd_gauss(r,[4.5 0.6]);

info = exp_5pdeer(t);
parIn = [info.parameters.default];

kappa = 0.4;
Bmodel = @(t) bg_exp(t,kappa);
K = exp_5pdeer(t,r,parIn,Bmodel);

V = K*P + whitegaussnoise(t,0.01);

[Vfit,Pfit] = fitsignal(V,t,r,'P',@bg_exp,@exp_5pdeer);

% Pass 1: signal is well fitted
pass(1) = all(abs(Vfit - V) < 3e-2);
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
   plot(t,V,'k.',t,Vfit,'r')
   legend('data','fit')
   xlabel('t [\mus]')
   ylabel('V(t)')
   grid on, axis tight, box on
end

end