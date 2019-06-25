function make_bas_Tikh()
%
% Computes kernel data for the DEER fit program in new format 2015
%
% G. Jeschke, 2015
%
dimension = 512;
ny0=52.04; % dipolar frequency at 1 nm for g=ge
dt = 0.008;
rmin = (4*dt*ny0/0.85)^(1/3);     % minimum detectable distance, 
                                  % the largest dipolar frequency
                                  % (shoulders of the Pake doublet is 85%
                                  % of the Nyquist frequency)
rmax = 7.5*(dimension*dt/2)^(1/3);  % maximum detectable distance, 
                                    % 2 microseconds dipolar evolution
                                    % time correspond to 6 nm, 25% margin
w0=2*pi*ny0; % angular frequencies
t=0:dt:(dimension-1)*dt; % time axis in µs
n=length(t); % length of time axis
r=linspace(rmin,rmax,dimension); % distance axis in nm
m=length(r); % length of distance axis
wd=zeros(1,m); % vector dipolar frequencies as a function of distance
kernel=zeros(m,n); % data array for kernel
fprintf(1,'Computing Tikhonov kernel');
fprintf(1,'%i data points in time domain\n',n);
fprintf(1,'%i data points in distance domain\n',m);
tic, % initialize computation time clock
for k=1:m
   rk=r(k); % current distance
   wdd=w0/rk^3; % current dipolar angular frequency
   wd(k)=wdd/(2*pi); % current dipolar frequency
   if mod(k,10)==0
        fprintf(1,'%5.1f%% of fit table done.\n',100*k/m);
   end
   for l=0:1000 % 1000 averages over cos(theta) angle (powder average)
      x=l/1000; % current theta angle
      ww=wdd*(3*x^2-1); % dipolar frequency at current theta angle
      kernel(k,:)=kernel(k,:)+cos(ww*t); % add kernel contribution
   end
end
toc, % output of computation time required
for k=1:m % loop form kernel normalization
   kernel(k,:)=kernel(k,:)/kernel(k,1); % normalize dipolar time evolution traces
end
L = get_l(length(r),2); % differential operator matrix for second derivative
[U,sm,X,V] = cgsvd(kernel',L);
save('pake_base_tikh_512.mat','kernel','r','t','U','sm','X','V','L'); % save kernel data
