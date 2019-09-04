function [err,data,maxerr] = test(opt,olddata)

%Test if selectmethod can identify that the optimal method is a two
%gaussian model as given as the input signal

Dimension = 200;
dt = 0.016;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [3 0.3 5 0.3 0.5];
P = rd_twogaussian(r,InputParam);

K = dipolarkernel(t,r);
S = K*P;

Models = {@rd_onegaussian,@rd_twogaussian,@rd_threegaussian};

%nm
[optimum1,metric1] = selectmodel(Models,S,r,K,'aicc');

%A
r = 10*time2dist(t);
[optimum2,metric2] = selectmodel(Models,S,r,K,'aicc');



err = optimum1~=optimum2;
data = [];
maxerr = [];


if opt.Display
figure(8),clf
hold on
plot(metric1)
plot(metric2)
end

end

