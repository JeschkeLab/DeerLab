function [pass,maxerr] = test(opt)

% Check that fitmultigauss() and fitparamodel() return the same answer

t = linspace(0,3.2,200);
r = time2dist(t);
InputParam = [3 0.5 4 0.5 0.4];
P = dd_gauss2(r,InputParam);
K = dipolarkernel(t,r);
S = K*P;
par0 = [2 0.1 5 0.1 0.5];
[~,P_FP] = fitparamodel(S,@dd_gauss2,r,K,par0);
P_MG = fitmultigauss(S,K,r,6);

% Pass: fitparamodel and fitmultigauss find the same solution
pass = all(abs(P_FP - P_MG) < 1e-9);

maxerr = max(abs(P_FP - P_MG));
 
if opt.Display
   plot(r,P,'k',r,P_MG,r,P_FP)
   legend('truth','multigauss','fitparamodel')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end