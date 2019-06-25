ra=2.0;
re=4.0;
noise=0.01;
dens=0.2;
dim=3;
depth=0.3;
fname='boxcar_20_40';

r=linspace(1.5,10,1000);
t=0:8:2400;
distr=zeros(size(r));
for k=1:length(distr),
    if r(k)>=ra & r(k)<=re,
        distr(k)=1;
    end;
end;
distr=0.01*distr/sum(distr);

deer=make_test_data(fname,r,distr,t,noise,dens,dim,depth);

figure(1); clf;
plot(r,distr,'b');
axis([1.5,8,-0.1*max(distr),1.1*max(distr)]);
axis off

figure(2); clf;
plot(t,deer,'k');