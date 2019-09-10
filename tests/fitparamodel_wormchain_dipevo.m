function [err,data,maxerr] = test(opt,oldata)


Dimension = 500;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [4 4.5];
L=InputParam(1);
Lp=InputParam(2);
kappa=Lp/L;
rn = r/L;
P=zeros(size(r));
%How Gunnar calculated the model
terms=2;
for pp=1:length(rn)
    G=0;
    crit=kappa*(1-rn(pp));
    if crit>0.2
        fac=2*kappa/(4*pi);
        for k=1:terms
            G=G+fac*pi^2*k^2*(-1)^(k+1)*exp(-kappa*pi^2*k^2*(1-rn(pp)));
        end
    elseif crit>0
        fac=kappa/(4*pi*2*sqrt(pi));
        for l=1:terms
            harg=(l-1/2)/sqrt(kappa*(1-rn(pp)));
            h2=4*harg^2-2;
            G=G+fac*1/(kappa*(1-rn(pp)))^(3/2)*exp(-(l-1/2)^2/(kappa*(1-rn(pp))))*h2;
        end
    end
    P(pp) = G;
end

P = P.';
P = P/sum(P)/mean(diff(r));

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

[FitParam,FitP] = fitparamodel(DipEvoFcn,@rd_wormchain,r,K,'Solver','lsqnonlin');
err(1) = any(abs(FitP - P)>5e-3);
err(2) = any(abs(FitParam - InputParam)>1e-2);
err = any(err);

maxerr = max(abs(FitP - P));
data = [];

if opt.Display
   figure(1),clf,
   subplot(121),hold on
   plot(t,DipEvoFcn,'b')
   plot(t,K*FitP,'r')
   subplot(122),hold on
   plot(r,P,'b')
   plot(r,FitP,'r')
end

end