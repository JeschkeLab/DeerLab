function [pass,maxerr] = test(opt)

%===============================================================================
% Make sure fitregparam handles custom regularization functional
%===============================================================================

N = 100;
dt = 0.008;
t = linspace(0,dt*N,N);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);
K = dipolarkernel(t,r);
V = K*P;

% Set optimal regularization parameter
alpha = 0.2;
L = regoperator(N,3);
RegFunctional = @(P)(1/2*norm(K*P-V)^2 + alpha^2/2*norm(L*P)^2);
Pfit = fitregmodel(V,K,r,RegFunctional,alpha,'Solver','fmincon');

deltaP = abs(Pfit-P);
pass = all(deltaP>0.1);
maxerr = max(deltaP);

 

if opt.Display
   	figure(8),clf
    hold on
    plot(r,P,'k') 
    plot(r,Pfit,'r')
    legend('model','fit');
end

end