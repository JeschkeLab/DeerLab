function [pass,maxerr] = test(opt)

%Test that using the multistart option helps finding the global minimum instead of the local one

t = linspace(0,5,300);
r1 = linspace(2,6,200);
P = rd_threegaussian(r1,[3 0.3 4 0.3 5 0.3 0.3 0.3]);
K = dipolarkernel(t,r1);
rng(5)
S = K*P + whitegaussnoise(t,0.01);
[~,Plocal] = fitparamodel(S,@rd_threegaussian,r1,K,'tolfun',1e-4);
[~,Pmulti] = fitparamodel(S,@rd_threegaussian,r1,K,'tolfun',1e-4,'multistart',50);

%Pass 1: solution with multi-start is better
pass(1) = max(abs(P - Pmulti)) < max(abs(P - Plocal));
%Pass 2: solution with multi-start fits the truth
pass(2) = all(abs(P - Pmulti) < 1e-1);

pass = all(pass);

maxerr = max(abs(P - Pmulti));
 

if opt.Display
plot(r1,P,'k',r1,Plocal,'r',r1,Pmulti,'b','linewidth',1)
grid on, axis tight
legend('truth','local minimum','global minimum')
xlabel('r [nm]')
ylabel('P(r) [nm^{-1}]')
end

end
