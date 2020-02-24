function [err,data,maxerr] = test(opt,olddata)


%Test that using the multistart option helps finding the global minimum
%instead of the local one

t = linspace(0,5,300);
r1 = linspace(2,6,200);
InputParam = [3 0.3 4 0.3 5 0.3 0.3 0.3];
P = rd_threegaussian(r1,InputParam);
K = dipolarkernel(t,r1);
rng(5)
S = K*P + whitegaussnoise(t,0.01);
[~,Pfit1] = fitparamodel(S,@rd_threegaussian,r1,K,'tolfun',1e-4);
r2 = linspace(2,6,200);
K = dipolarkernel(t,r2);
[~,Pfit2] = fitparamodel(S,@rd_threegaussian,r2,K,'tolfun',1e-4,'multistart',50);

err(1) = max(abs(P - Pfit2)) > max(abs(P - Pfit1)>1e-5);
err(2) = any(abs(P - Pfit2) > 1e-1);

err = any(err);
maxerr = max(abs(P - Pfit2));
data = [];

if opt.Display
clf
plot(r1,P,'k',r1,Pfit1,'r',r2,Pfit2,'b','linewidth',1)
grid on, axis tight
legend('truth','Local','Global')
xlabel('r [nm]')
ylabel('P(r) [nm^{-1}]')
set(gca,'fontsize',13)
end

end
