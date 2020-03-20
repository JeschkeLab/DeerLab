function [pass,maxerr] = test(opt)

% Check error control of nosielevel() towards wrong inputs

rng(1)
t = linspace(0,3,200);
r = linspace(2,6,100);
P = dd_onegaussian(r,[3,0.3]);
K = dipolarkernel(t,r);
S = dipolarsignal(t,r,P,'noiselevel',0.05,'phase',pi/2);

alpha = 0.5;

% Pass 1: passing a signal with complex values
try
    obir(S,K,r,'tikh',alpha);
    pass = false;
catch
    pass = true;
end

maxerr = NaN;
 

end