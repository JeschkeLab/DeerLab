function [err,data,maxerr] = test(opt,olddata)

%Test if selectmethod can identify that the optimal method is a two
%rice model as given as the input signal
Dimension = 300;
dt = 0.016;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [3 0.2 5.5 0.7 0.5];
P = rd_tworice(r,InputParam);
P = P/sum(P)/mean(diff(r));

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

Models = {@rd_onerice,@rd_tworice,@rd_threerice};

[optimum,metric] = selectmodel(Models,DipEvoFcn,r,K,'aicc','Solver','lsqnonlin');

err = optimum~=2;
data = [];
maxerr = 0;


if opt.Display
figure(8),clf
plot(metric)
end

end

