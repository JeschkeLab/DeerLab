ra=2.0;
re=4.0;
noise=0.01;
dens=0.2;
dim=3;
depth=0.3;
fname='sawtooth_20_40';

r=linspace(1.5,10,1000);
[mo,pa]=min(abs(r-ra*ones(size(r))));
[mo,pe]=min(abs(r-re*ones(size(r))));
rm=(ra+re)/2;
[mo,pm]=min(abs(r-rm*ones(size(r))));
t=0:8:2400;
distr=zeros(size(r));
for k=pa:pm,
   distr(k)=(k-pa)/(pm-pa);
end;
for k=pm+1:pe,
   distr(k)=1-(k-pm)/(pe-pm);
end;
distr=0.01*distr/sum(distr);

deer=make_test_data(fname,r,distr,t,noise,dens,dim,depth);

figure(1); clf;
plot(r,distr,'b');
axis([1.5,8,-0.1*max(distr),1.1*max(distr)]);
axis off

figure(2); clf;
plot(t,deer,'k');