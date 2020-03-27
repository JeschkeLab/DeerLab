function [pass,maxerr] = test(opt)

% Check that fitbackground() works with negative times

t = linspace(-0.5,5,300);
r = time2dist(t);
P = dd_onegauss(r,[4 0.4]);
K = dipolarkernel(t,r);
rng(2)
S = K*P + whitegaussnoise(300,0.02);
lam = 0.2;
k = 0.3;
B = exp(-k*(abs(t)))';
F = (1-lam) + lam*S;
V = F.*B;
B = (1-lam).*B;
tstart = backgroundstart(V,t,@bg_exp);
[Bfit,lambdafit] = fitbackground(V,t,@bg_exp,tstart);
Bfit = (1-lambdafit).*Bfit;

error = abs(B - Bfit);

% Pass: backgrund is well fitted
pass = all(error < 1e-2);

maxerr = max(error);
 

if opt.Display
    plot(t,B,t,Bfit)
    legend('truth','fit')
    xlabel('t [\mus]')
    ylabel('B(t)')
    grid on, axis tight, box on
end

end