ra=3.0;
re=5.0;
noise=0.01;
dens=0.2;
dim=3;
depth=0.3;
fname='halftooth_30_50_noise';

r=linspace(1.5,10,1000);
[mo,pa]=min(abs(r-ra*ones(size(r))));
[mo,pe]=min(abs(r-re*ones(size(r))));
t=0:8:2400;
distr=zeros(size(r));
for k=pa:pe,
   distr(k)=(k-pa)/(pe-pa);
end;
distr=0.01*distr/sum(distr);

deer=make_test_data(fname,r,distr,t,noise,dens,dim,depth);

figure(1); clf;
plot(r,distr,'b');
axis([1.5,8,-0.1*max(distr),1.1*max(distr)]);
axis off

figure(2); clf;
plot(t,deer,'k');