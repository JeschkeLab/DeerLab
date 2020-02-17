function [err,data,maxerr] = test(opt,olddata)

%===============================================================================
% Check whether kernel is constructed properly for negative times
%===============================================================================

dt = 0.020;
n = 50;
t = (-n:n)*dt;

r = 3;
K = dipolarkernel(t,r,'Method','grid');

negK = K(1:n+1);
posK = K(n+1:end);

delta = abs(negK - posK(end:-1:1));

err(1) = any(delta(:)>1e-12);
err = any(err);
maxerr = max(delta(:));
data = [];

if opt.Display
  plot(t,K,t,K(end:-1:1));
  legend('as calculated','flipped');
end

end