function deer=make_test_data(fname,r0,distr0,t0,noise,dens,dim,depth)
%
% Creates a DEER test data set with file name fname[.dat] in ASCII format
% (1st column time in ns, 2nd column intensity)
% from a given distance distribution for a given time axis provided in
% variable t (in ns)
% r0        distance axis
% distr0    distribution
% t0        time axis
% noise     noise level (signal is normalized at maximum)
% dens      background density
% dim       background dimension
% depth     modulation depth

t0=t0/1000; % conversion to microseconds

% load kernel
load('pake_base40.mat');
kernel=base-ones(size(base)); % kernel format for pcf2deer

% interpolate distribution to kernel distance axis and normalize
distr=get_std_distr(r0,distr0,r);
distr=0.01*distr/sum(distr);

% transform distribution to form factor
deer0=pcf2deer(distr,kernel,r,t);
deer0=deer0-0.99*ones(size(deer0));

% interpolate to input time axis
deer1=interp1(t,deer0,t0);

% normalize
deer1=deer1/max(deer1);

% figure(13); clf;
% plot(t0,deer1,'k');

% add unmodulated part
deer1=(1-depth)*ones(size(deer1))+depth*deer1;
% hold on;
% plot(t0,deer1,'r');

% make background function B(t)
targ=t0.^(dim/3);
bckg=exp(-dens*targ);

% convolute with background
deer=bckg.*deer1;

% add (pseudorandom) noise
deer=deer+noise*randn(size(deer));

data=[1000*t0' deer'];
save([fname '.dat'],'data','-ascii');

data2=[r0' distr0'];
save([fname '_theor.distr'],'data2','-ascii');

