function [pass,maxerr] = test(opt)

%Test if selectmethod can identify that the optimal method is a two
%gaussian model as given as the input signal

Dimension = 200;
dt = 0.016;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [3 0.3 5 0.3 0.5];
P = rd_twogaussian(r,InputParam);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

Models = {@rd_onegaussian,@rd_twogaussian,@rd_threegaussian};

[optimum1,metric1] = selectmodel(Models,DipEvoFcn,r,K,'aicc');
[optimum2,metric2] = selectmodel(Models,DipEvoFcn,r,K,'aic');
[optimum3,metric3] = selectmodel(Models,DipEvoFcn,r,K,'bic');
[optimum4,metric4] = selectmodel(Models,DipEvoFcn,r,K,'rmsd');

err(1) = optimum1~=optimum2;
err(2) = optimum2~=optimum3;
pass = all(err);
 
maxerr = NaN;


if opt.Display
figure(8),clf
subplot(221)
plot(metric1)
subplot(222)
plot(metric2)
subplot(223)
plot(metric3)
subplot(224)
plot(metric4)

end

end

