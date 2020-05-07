function [pass,maxerr] = test(opt)

% Check that fitregmodel handles custom regularization functional

t = linspace(0,3.2,200);
r = linspace(2,6,100);
P = dd_gauss(r,[3,0.5]);
K = dipolarkernel(t,r);
S = K*P;
alpha = 0.05;
L = regoperator(r,3);
RegFunctional = @(P)(1/2*norm(K*P - S)^2 + alpha^2/2*norm(L*P)^2);
Pfit = fitregmodel(S,K,r,RegFunctional,alpha,'Solver','fmincon');
deltaP = abs(Pfit - P);

% Pass: distribution is well fitted
pass = all(deltaP < 1e-2);

maxerr = max(deltaP);

if opt.Display
   plot(r,P,r,Pfit)
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end