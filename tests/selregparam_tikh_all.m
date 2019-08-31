function [err,data] = test(opt,olddata)

%=======================================
% Check regparamrange.m
%=======================================

Dimension = 100;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
Distribution = rd_onegaussian(r,[3,0.5]);
Distribution = Distribution/sum(Distribution);

K = dipolarkernel(t,r);
RegMatrix = regoperator(Dimension,2);
DipEvoFcn = K*Distribution;

RegParamSet = regparamrange(K,RegMatrix);
OptParam2 = selregparam(RegParamSet,DipEvoFcn,K,RegMatrix,'tikhonov','all','NonNegConstrained',false,'NoiseLevel',0.05);

%Accept testif all values are the same (should be as there is no noise)
err = length(OptParam2)~=15;
data = [];



end