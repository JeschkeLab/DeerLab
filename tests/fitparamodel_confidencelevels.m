function [pass,maxerr] = test(opt)

% Test that confidence levels can be adjusted via option

t = linspace(0,8,300);
r = linspace(1,8,300);
parIn  = [4.5,0.5];
P = dd_onegauss(r,parIn);

K = dipolarkernel(t,r);
S = K*P + whitegaussnoise(t,0.01);

par0 = [2 0.2];
[parFit1,Pfit,parCI1] = fitparamodel(S,@dd_onegauss,r,K,par0,'confidencelevel',0.95);
[parFit2,~,parCI2] = fitparamodel(S,@dd_onegauss,r,K,par0,'confidencelevel',0.50);
[parFit3,~,parCI3] = fitparamodel(S,@dd_onegauss,r,K,par0,'confidencelevel',0.25);


% Pass 1-2: confidence intervals behave as expected
pass(1) = all(all(abs(parFit1 - parCI1.') > abs(parFit2 - parCI2.')));
pass(2) = all(all(abs(parFit2 - parCI2.') > abs(parFit3 - parCI3.')));

maxerr = NaN;

 
if opt.Display
   plot(r,P,'k',r,Pfit,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end