function sim=sim_dir_Tikh(distr)
%
% Computes kernel data for the DEER fit program
%
% G. Jeschke, 2002
%
dimension = 1024;
rmin = 1.5;
rmax = 15;
ny0=52.04; % dipolar frequency at 1 nm for g=ge
w0=2*pi*ny0; % angular frequencies
dt = 0.008;
t=0:dt:(dimension-1)*dt; % time axis in µs
n=length(t); % length of time axis
r=linspace(rmin,rmax,dimension); % distance axis in nm
m=length(r); % length of distance axis
wd=zeros(1,m); % vector dipolar frequencies as a function of distance
fprintf(1,'%i data points in time domain\n',n);
fprintf(1,'%i data points in distance domain\n',m);
sim = zeros(1,n);
tic, % initialize computation time clock
for k=1:m,
   rk=r(k); % current distance
   wdd=w0/rk^3; % current dipolar angular frequency
   wd(k)=wdd/(2*pi); % current dipolar frequency
   fprintf(1,'%5.1f%% of fit table done.\n',100*k/m);
   for l=0:1000, % 1000 averages over cos(theta) angle (powder average)
      x=l/1000; % current theta angle
      ww=wdd*(3*x^2-1); % dipolar frequency at current theta angle
      sim = sim + distr(k)*cos(ww*t); % add kernel contribution
   end;
end;
toc, % output of computation time required

