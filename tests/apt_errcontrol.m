function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Test data preparation function
%======================================================

%Parameters
Dimension = 200;
t = linspace(-0.6,2,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

aptK = aptkernel(t);
try
    apt(DipEvoFcn,aptK,[0.05 125]);
    err = true;
catch
    err = false;
end

K = rand(100,100);
try
    apt(DipEvoFcn,K,2);
    err = true;
catch
    err = false;
end

S = DipEvoFcn + 1i*DipEvoFcn;
try
    apt(S,aptK,2);
    err = true;
catch
    err = false;
end

err = any(err);
data = [];
maxerr = NaN;

end