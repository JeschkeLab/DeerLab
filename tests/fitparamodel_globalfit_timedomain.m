function [err,data,maxerr] = test(opt,olddata)

Ntime1 = 100;
Ndist = 200;
rng(2)
dt = 0.008;
t1 = linspace(0,dt*Ntime1,Ntime1);
[~,rmin,rmax] = time2dist(t1);
r = linspace(rmin,rmax,Ndist);

P = rd_twogaussian(r,[2,0.3,4,0.3,0.5]);

K1 = dipolarkernel(t1,r);
S1 = K1*P;
noise = whitegaussnoise(length(S1),0.03);
S1 = S1 + noise;

Ntime2 = 200;
t2 = linspace(0,dt*Ntime2,Ntime2);
K2 = dipolarkernel(t2,r);
S2 = K2*P;
noise = whitegaussnoise(length(S2),0.05);
S2 = S2 + noise;

Ntime3 = 300;
t3 = linspace(0,dt*Ntime3,Ntime3);
K3 = dipolarkernel(t3,r);
S3 = K3*P;
noise = whitegaussnoise(length(S3),0.1);
S3 = S3 + noise;
rng('default')

Ss = {S1,S2,S3};

tmodel = t3;
info = rd_twogaussian();
mymodel = @(t,param)dipolarkernel(t,r)*rd_twogaussian(r,param);


InitialGuess = [2 0.1 5 0.1 0.5];
range = [info.parameters(:).range];
upper = range(2:2:end);
lower = range(1:2:end-1);

paramglobal = fitparamodel(Ss,mymodel,{t1,t2,t3},InitialGuess,'Upper',upper,'lower',lower);

param3 = fitparamodel(S3,mymodel,tmodel,InitialGuess,'Upper',upper,'lower',lower);

Pfit = rd_twogaussian(r,paramglobal);
Dist3 = rd_twogaussian(r,param3);


normResult = norm(P - Pfit);
norm3 = norm(P - Dist3);

err(1) = any(normResult > [norm3]);
err = any(err);
data = [];
maxerr = max(P - Pfit);

if opt.Display
    figure(8),clf
    subplot(121)
    hold on
    plot(t1,S1,'k')
    plot(t1,K1*Pfit,'r')
    plot(t2,S2+1,'k')
    plot(t2,K2*Pfit + 1,'r')
    plot(t3,S3+2,'k')
    plot(t3,K3*Pfit + 2,'r')
    subplot(122)
    hold on
    plot(r,P,'k')
    plot(r,Pfit,'r')
    plot(r,Dist3,'b--')
end

end