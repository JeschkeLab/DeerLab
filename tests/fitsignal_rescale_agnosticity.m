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

[Vfit1,Pfit1] = fitsignal(V*scale,t,r,@dd_gauss,@bg_exp,@ex_4pdeer);
[Vfit2,Pfit2] = fitsignal(V,t,r,@dd_gauss,@bg_exp,@ex_4pdeer);

%Pass 1-2: distance distributions are well fitted
pass(1) = all(abs(Pfit1 - P) < 2e-1);
pass(2) = all(abs(Pfit2 - Pfit1) < 1e-2);

%Pass 3-4: dipolar signals are well fitted
pass(3) = all(abs(Vfit1 - V*scale) < scale*2e-2);
pass(4) = all(abs(Vfit2 - V) < 2e-2);

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