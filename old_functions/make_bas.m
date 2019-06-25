function make_bas()
%
% Computes kernel data for the DEER fit program
%
% G. Jeschke, 2002,2019
%
fprintf(1,'Computing simulation kernel. Please wait.\n');
ny0=52.04; % dipolar frequency at 1 nm for g=ge
w0=2*pi*ny0; % angular frequencies
t=0:0.008:16.384-0.008; % time axis in µs
n=length(t); % length of time axis
r=1.5:0.02:40; % distance axis in nm
m=length(r); % length of distance axis
wd=zeros(1,m); % vector dipolar frequencies as a function of distance
base=zeros(m,n); % data array for kernel
tic, % initialize computation time clock
for k=1:m
   rk=r(k); % current distance
   wdd=w0/rk^3; % current dipolar angular frequency
   wd(k)=wdd/(2*pi); % current dipolar frequency
   if mod(k,10) == 0
    fprintf(1,'%5.1f of fit table done\n',100*k/m);
   end
   for x=0:0.001:1 % 1000 averages over theta angle (powder average)
      ww=wdd*(3*x^2-1); % dipolar frequency at current theta angle
      base(k,:)=base(k,:)+cos(ww*t); % add kernel contribution
   end
end
toc, % output of computation time required
for k=1:m % loop form kernel normalization
   base(k,:)=base(k,:)/base(k,1); % normalize dipolar time evolution traces
end
save('pake_base40.mat','base','r','t','wd'); % save kernel data
fprintf(1,'Simulation kernel computed.\n');

