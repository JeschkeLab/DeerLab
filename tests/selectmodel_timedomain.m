function [err,data,maxerr] = test(opt,olddata)

%Test if selectmethod can identify that the optimal method is a two
%gaussian model as given as the input signal

Dimension = 200;
dt = 0.016;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);

K = dipolarkernel(t,r);
B = td_strexp(t,[0.2 4]);

Models = {@td_exp,@td_strexp,@td_prodstrexp};

[optimum1,metric] = selectmodel(Models,B,t,'aicc');
optimum2 = selectmodel(Models,B,t,'aic');
optimum3 = selectmodel(Models,B,t,'bic');

err(1) = optimum1~=optimum2;
err(2) = optimum2~=optimum3;
err = any(err);
data = [];
maxerr = [];


if opt.Display
figure(8),clf
plot(metric)
end

end

