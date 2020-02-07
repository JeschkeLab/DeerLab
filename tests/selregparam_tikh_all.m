function [err,data,maxerr] = test(opt,olddata)

Dimension = 100;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);
P = P/sum(P);

K = dipolarkernel(t,r);
S = K*P;

OptParam1 = selregparam(S,K,r,'tikhonov','all','NonNegConstrained',false,'NoiseLevel',0.05);
OptParam2 = selregparam(S,K,r,'tikhonov',{'all'},'NonNegConstrained',false,'NoiseLevel',0.05);

%Accept testif all values are the same (should be as there is no noise)
err(1) = length(OptParam2)~=14;
err(2) = length(OptParam1)~=14;
maxerr = NaN;
data = [];



end