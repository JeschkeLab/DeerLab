function [pass,maxerr] = test(opt)

% Check that fitsignal() does pre-processing if required

rng(1)

t = linspace(0,5,100);
r = linspace(2,6,150);
P = dd_gauss(r,[4.5 0.6]);

info = ex_4pdeer(t);
parIn = [info.parameters.default];
pathinfo = ex_4pdeer(t,parIn);

kappa = 0.4;
Bmodel = @(t,lam) bg_exp(t,kappa,lam);
V = dipolarsignal(t,r,P,pathinfo,Bmodel,'noiselevel',0.01,'scale',1e5,'phase',pi/6);

[Vfit,Pfit] = fitsignal(V,t,r,'P',@bg_exp,@ex_4pdeer);


V = correctphase(V);
V = correctscale(V,t,0.5);

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