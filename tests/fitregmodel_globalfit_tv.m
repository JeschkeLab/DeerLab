function [err,data,maxerr] = test(opt,olddata)

warning('off','all')

Ntime1 = 100;
Ndist = 200;

dt = 0.008;
t1 = linspace(0,dt*Ntime1,Ntime1);
[~,rmin,rmax] = time2dist(t1);
r = linspace(rmin,rmax,Ndist);

Distribution = rd_onegaussian(r,[2,0.3]) + rd_onegaussian(r,[4,0.3]);
Distribution = Distribution/sum(Distribution)/mean(diff(r));
K1 = dipolarkernel(t1,r);
S1 = K1*Distribution;
noise = whitenoise(length(S1),0.03);
S1 = S1 + noise;

Ntime2 = 200;
t2 = linspace(0,dt*Ntime2,Ntime2);
K2 = dipolarkernel(t2,r);
S2 = K2*Distribution;
noise = whitenoise(length(S2),0.05);
S2 = S2 + noise;

Ntime3 = 300;
t3 = linspace(0,dt*Ntime3,Ntime3);
K3 = dipolarkernel(t3,r);
S3 = K3*Distribution;
noise = whitenoise(length(S3),0.1);
S3 = S3 + noise;


L = regoperator(Ndist,2);
%Set optimal regularization parameter (found numerically lambda=0.13)
regparam = 1;

Ss = {S1,S2,S3};
Ks = {K1,K2,K3};

Result = fitregmodel(Ss,Ks,r,L,'tv',regparam,'Solver','fnnls');
Dist1 = fitregmodel(S1,K1,r,L,'tv',regparam,'Solver','fnnls');
Dist2 = fitregmodel(S2,K2,r,L,'tv',regparam,'Solver','fnnls');
Dist3 = fitregmodel(S3,K3,r,L,'tv',regparam,'Solver','fnnls');

normResult = norm(Distribution - Result);
norm1 = norm(Distribution - Dist1);
norm2 = norm(Distribution - Dist2);
norm3 = norm(Distribution - Dist3);

err(1) = any(normResult > [norm2 norm3]);
err = any(err);
data = [];
maxerr = normResult;

if opt.Display
figure(8),clf
subplot(121)
hold on
plot(t1,S1,'k')
plot(t2,S2+1,'k')
plot(t3,S3+2,'k')
plot(t1,K1*Result,'r')
plot(t2,K2*Result + 1,'r')
plot(t3,K3*Result + 2,'r')
subplot(122)
hold on
plot(r,Distribution,'k')
plot(r,Result,'r')
plot(r,Dist1,'g--')
plot(r,Dist2,'b--')
plot(r,Dist3,'r--')
end

end