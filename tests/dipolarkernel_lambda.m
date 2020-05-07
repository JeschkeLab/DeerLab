function [pass,maxerr] = test(opt)

% Check that dipolar kernel with modulation depth works

t = linspace(0,4,150); % us
r = time2dist(t); % nm
P = dd_gauss(r,[3,0.3]);
K = dipolarkernel(t,r);

lam = 0.4;
Fref = (1-lam) + lam*K*P;
Klam = dipolarkernel(t,r,lam);
Ffit = Klam*P;

err = abs(Fref-Ffit);
maxerr = max(err);
pass = maxerr < 1e-10;
 
% Plot if requested
if opt.Display
   plot(t,Fref,'k',t,Ffit,'r')
   legend('truth','K*P')
   xlabel('t [\mus]')
   ylabel('V(t)')
   grid on, axis tight, box on
end


end