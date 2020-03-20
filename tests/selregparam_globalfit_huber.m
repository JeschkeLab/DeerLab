function [pass,maxerr] = test(opt)

rng(1)

t1 = linspace(0,4,100);
t2 = linspace(0,3,200);
t3 = linspace(0,4,300);

r = linspace(2,5,200);
P = dd_twogaussian(r,[3,0.4,3.5,0.4,0.3]);

K1 = dipolarkernel(t1,r);
S1 = dipolarsignal(t1,r,P,'noiselevel',0.03);
K2 = dipolarkernel(t2,r);
S2 = dipolarsignal(t2,r,P,'noiselevel',0.02);
K3 = dipolarkernel(t3,r);
S3 = dipolarsignal(t3,r,P,'noiselevel',0.02);

Ss = {S1,S2,S3};
Ks = {K1,K2,K3};

alpha = selregparam(Ss,Ks,r,'huber','aic');

Pfit = fitregmodel(Ss,Ks,r,'huber',alpha);

% Pass: the distribution is well fitted with the optimized parameter
pass = all(abs(Pfit - P) < 0.5);
 
maxerr = max(abs(Pfit - P));

if opt.Display
   plot(r,P,r,Pfit)
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end