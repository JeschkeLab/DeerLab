function [pass,maxerr] = test(opt)

% Check that fitregmodel() works with non-square kernel matrices

t = linspace(0,2,200);
r = linspace(2,6,100);
P = rd_onegaussian(r,[3,0.5]);
K = dipolarkernel(t,r);
S = K*P;

Pfit = fitregmodel(S,K,r,'tikhonov','aic','Solver','fnnls');

% Pass 1: the distribution is well fitted
pass(1) = all(abs(Pfit - P) < 1e-3);
% Pass 2: the distance-dimension has the right size
pass(2) = length(Pfit) == numel(r);
% Pass 3: the time-dimension has the right size
pass(3) = length(K*Pfit) == numel(t);
pass  = any(pass);

maxerr = max(abs(Pfit - P));
 
if opt.Display
   plot(r,P,r,Pfit)
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end