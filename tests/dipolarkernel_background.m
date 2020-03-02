function [pass,maxerr] = test(opt)

% Test kernel construction with modulation depth and background

t = linspace(0,3,80);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);
B = td_exp(t,0.5);
K = dipolarkernel(t,r);

S  = K*P;
lam = 0.25;
V = (1 - lam + lam*S).*B;

KB = dipolarkernel(t,r,lam,B);
Vfit  = KB*P;

% Pass: the kernel transform the distribution into correct signal
pass = all(abs(Vfit - V) < 1e-10);

maxerr = max(abs(Vfit - V));
 

if opt.Display
   plot(t,V,'k',t,Vfit,'r',t,(1-lam)*B,'r--')
   legend('truth','B','K_B*P')
   xlabel('t [\mus]')
   ylabel('V(t)')
   grid on, axis tight, box on
end

end