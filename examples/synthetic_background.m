% you need once, with DeerAnalysis on the path
%load('pake_base40.mat');
%kernel=base-ones(size(base)); % kernel format for pcf2deer
distr=r.^2;
distr=50*distr/sum(distr);
deer=pcf2deer(distr,kernel,r,t);
deer=deer+0.005*randn(1,1024); % some pseudo-noise

data=[1000*t(1:512)',deer(1:512)'];
save synthetic_3D_background.dat data -ascii

figure(1); clf;
plot(t,deer,'k');
hold on

distr=r;
distr=15*distr/sum(distr);
deer=pcf2deer(distr,kernel,r,t);
deer=deer+0.005*randn(1,1024); % some pseudo-noise

data=[1000*t(1:512)',deer(1:512)'];
save synthetic_2D_background.dat data -ascii

plot(t,deer,'b');


distr=ones(size(r));

distr=5*distr/sum(distr);
deer=pcf2deer(distr,kernel,r,t);
deer=deer+0.005*randn(1,1024); % some pseudo-noise

data=[1000*t(1:512)',deer(1:512)'];
save synthetic_1D_background.dat data -ascii

plot(t,deer,'r');

