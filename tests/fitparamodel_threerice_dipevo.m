function [err,data,maxerr] = test(opt,oldata)


Dimension = 300;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [2 0.4 3.5 0.3 5 0.3 0.3 0.3];
P = rd_threerice(r,InputParam);
P = P/sum(P)/mean(diff(r));

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

[~,FitP] = fitparamodel(DipEvoFcn,@rd_threerice,r,K,0.75*InputParam,'solver','lsqnonlin');
err = any(abs(FitP - P)>1e-2);

maxerr = max(abs(FitP - P));
data = [];

if opt.Display
   figure(1),clf
   subplot(121)
   hold on
   plot(t,DipEvoFcn,'b')
   plot(t,K*FitP,'r')
   subplot(122)
   hold on
   plot(r,P,'b')
   plot(r,FitP,'r')
end

end