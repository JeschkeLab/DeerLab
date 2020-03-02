function [pass,maxerr] = test(opt)

% Check that dipolar kernel with modulation depth works

t = linspace(0,5,50); % us
r = time2dist(t); % nm
P = rd_onegaussian(r,[3,0.5]);
K = dipolarkernel(t,r);
lam = 0.4;
F = 1 -lam + lam*K*P;
K = dipolarkernel(t,r,lam);
Ffit = K*P;

% Pass: the kernel transform into the correct signal
pass = all(abs(F - Ffit) < 1e-10);

maxerr = max(abs(F - Ffit));
 
%Plot if requested
if opt.Display
   plot(t,V,'k',t,Ffit,'r')
   legend('truth','K*P')
   xlabel('t [\mus]')
   ylabel('V(t)')
   grid on, axis tight, box on
end


end