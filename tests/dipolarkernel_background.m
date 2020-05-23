function [pass,maxerr] = test(opt)

% Test kernel construction with modulation depth and background

t = linspace(0,3,80);
r = time2dist(t);
P = dd_gauss(r,[3,0.5]);
B = bg_exp(t,0.5);
K = dipolarkernel(t,r);

S  = K*P;
lam = 0.25;
V0 = (1 - lam + lam*S).*B;

KB = dipolarkernel(t,r,lam,B);
V  = KB*P;

maxerr = max(abs(V-V0));
pass = maxerr < 1e-5;

if opt.Display
   plot(t,V0,'k',t,V,'r',t,(1-lam)*B,'r--')
   legend('truth','B','K_B*P')
   xlabel('t [\mus]')
   ylabel('V(t)')
   grid on, axis tight, box on
end

end