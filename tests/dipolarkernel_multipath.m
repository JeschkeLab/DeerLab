function [pass,maxerr] = test(opt)

% Check that multi-pathway kernels are properly generated

r = linspace(2,6,50);

% Example timings for five-pulse DEER (all in us)
t1 = linspace(0,10,300);
t2 = 0.3;
tau1 = 4.24;
tau2 = 4.92;
t = (tau1 + tau2) - (t1 + t2);

% Pathway amplitudes and zero times
prob = 0.8;
lambda = [1-prob; prob^2; prob*(1-prob)];
T0 = [NaN; 0; tau2-t2];

K = dipolarkernel(t,r,[lambda T0]);

unmodulated = isnan(T0);
Kref = sum(lambda(unmodulated));
dr = mean(diff(r));
for p = 1:numel(lambda)
  if unmodulated(p), continue; end
  Kref = Kref + lambda(p)*dipolarkernel(t-T0(p),r)/dr;
end
Kref = Kref*dr;

maxerr = abs(max(K(:) - Kref(:)));

% Pass: the kernel returned by the function is equal to the reference
pass = maxerr < 1e-3;
 
if opt.Display
   plot(t,Kref(:,25),'k',t,K(:,25),'r')
   legend('reference','output')
   xlabel('t [\mus]')
   ylabel('K(t,4nm)')
   grid on, axis tight, box on
end

end
