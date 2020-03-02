function [pass,maxerr] = test(opt)

% Test the cosntruction of a mixed model of different basis functions

t = linspace(0,3,200);
r = linspace(2,6,100);
parIn1 = [2.5 0.3];
P1 = rd_onegaussian(r,parIn1);
parIn2 = [3.5 0.3];
P2 = rd_onegaussian(r,parIn2);
parIn3 = [4.5 0.3];
P3 = rd_onerice(r,parIn3);
P = 0.4*P2 + 0.3*P1 + 0.3*P3;

mixedModel = mixmodels({@rd_onegaussian,@rd_onegaussian,@rd_onerice});
parInMix = [0.3 0.4 parIn1 parIn2 parIn3];
Pmix = mixedModel(r,parInMix);

% Pass: the models have been mixed properly
pass = all(abs(Pmix - P) < 1e-8);

maxerr = max(abs(Pmix - P));
 
if opt.Display
   plot(r,P,'k',r,Pmix)
   legend('truth','mix')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end