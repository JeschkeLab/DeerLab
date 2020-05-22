function [pass,maxerr] = test(opt)

% Check that the different selregparam options work

t = linspace(0,3,200);
r = linspace(2,6,100);
P = dd_gauss(r,[3,0.5]);
K = dipolarkernel(t,r);
S = K*P;

alphaopt1 = selregparam(S,K,r,'tikhonov',{'aic','gml','gcv'},'NoiseLevel',1e-5);
alphaopt2 = selregparam(S,K,r,'tikhonov',{'aic','gml','gcv'},'RegOrder',2);
alphaopt3 = selregparam(S,K,r,'tikhonov',{'aic','gml','gcv'},'TolFun',1e-10);
alphaopt4 = selregparam(S,K,r,'huber',{'aic','gml','gcv'});
alphaopt5 = selregparam(S,K,r,'huber',{'aic','gml','gcv'},'HuberParameter',1.35);
alphaopt6 = selregparam(S,K,r,'huber',{'aic','gml','gcv'},'GlobalWeights',1);


% Pass 1-6: the different options wornd the results are similar
pass(1) = all(abs(alphaopt1 - alphaopt2) < 1e-2);
pass(2) = all(abs(alphaopt1 - alphaopt3) < 1e-2);
pass(3) = all(abs(alphaopt1 - alphaopt3) < 1e-2);
pass(4) = all(abs(alphaopt4 - alphaopt5) < 1e-2);
pass(5) = all(abs(alphaopt4 - alphaopt6) < 1e-2);

pass = all(pass);
 
maxerr = NaN;

end