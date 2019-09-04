function [err,data,maxerr] = test(opt,olddata)

Ntime1 = 100;
Ndist = 200;
dt = 0.008;
t1 = linspace(0,dt*Ntime1,Ntime1);
[~,rmin,rmax] = time2dist(t1);
r = linspace(rmin,rmax,Ndist);
P = rd_twogaussian(r,[2,0.3,3.5,0.3,0.5]);

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


L = regoperator(Ndist,2);

Ss = {S1,S2,S3};
Ks = {K1,K2,K3};


RegParamRange = logspace(-3,4,60);
[OptRegParam,fun] = selregparam(RegParamRange,Ss,Ks,L,'tikhonov',{'aic','aicc'});

Result1 = fitregmodel(Ss,Ks,r,L,'tikhonov',OptRegParam(1),'Solver','fnnls');
Result2 = fitregmodel(Ss,Ks,r,L,'tikhonov',OptRegParam(2),'Solver','fnnls');

Result = mean([Result1 Result2],2);
stdDist = std([Result1 Result2],1,2);

err = any(abs(Result-P) > 0.5);
data = [];
maxerr = max(abs(Result-P));

if opt.Display
figure(8),clf
subplot(121)
hold on
plot(RegParamRange,fun{1},'.b')
subplot(122)
hold on
plot(r,P,'k')
Result = Result';
stdDist = stdDist';
f = fill([r fliplr(r)],[Result+stdDist fliplr(Result-stdDist)],'b','LineStyle','none');
f.FaceAlpha = 0.5;
plot(r,Result,'b')
end

end