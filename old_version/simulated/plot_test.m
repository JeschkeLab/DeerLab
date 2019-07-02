fname='gaussian_20_5.dat';
data=load(fname);
t=data(:,1)/1000;
distr=data(:,2);
figure(4); clf;
plot(t,distr,'b');
