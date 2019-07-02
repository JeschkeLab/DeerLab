function [r,distr,sim]=APT(handles,tdip,dipevo)
%
% Pake transformation
%

base=handles.APT_kernel;
crosstalk=handles.APT_crosstalk;
trcnorm=handles.APT_norm;
t=handles.APT_t;
ny=handles.APT_ny;
r=handles.Pake_r;
kernel=handles.Pake_kernel;
fwdt=handles.DDS;

[m,n]=size(base); % size of kernel
dt=tdip(2)-tdip(1); % time increment
stretch=(dt/0.008)^(1/3);
nt=length(t); % length of time axis
spc2=zeros(1,m); % initialize distribution
td=zeros(1,n);
td(1:length(tdip))=dipevo;
tdx=td.*t; % eqn [21]
for k=1:m % sum in eqn [21]
  spc2(k)=spc2(k)+sum(base(k,:).*tdx)/trcnorm(k); 
end
spc3=crosstalk\spc2'; % crosstalk correction, eqn [22]

ny2r=ny; % initialize distance axis (mapping of dipolar frequencies to distances)
for k=1:length(ny2r),
    ny2r(k)=(52.04/ny2r(k))^(1/3);
end;
ny2r=stretch*ny2r;

ny2re=ny/2; % initialize distance axis (mapping of dipolar frequencies to distances)
for k=1:length(ny2re),
    ny2re(k)=(52.04/ny2re(k))^(1/3);
end;
ny2re=stretch*ny2re;

spc4=interp1(ny2r,spc3,ny2re,'pchip',0); % initialize data for distance domain smoothing
spc4=spc4';
spc5=spc4;
for k=1:length(spc4), % smoothing, convolute with...
    ract=ny2re(k);
    filt=(ny2re-ract*ones(1,length(ny2re)))/fwdt; % ... Gaussian line as filter
    filt=exp(-filt.*filt);
    filt=filt';
    spc5(k)=sum(filt.*spc4)/sum(filt); % normalization!
end;

for k=1:m,
    spc5(k)=spc5(k)/(ny2re(k))^4; % keep constant integral after mapping to distances
end;

spc5=spc5/max(spc5); % renormalize distribution


distr=get_std_distr(ny2re,spc5,r);
distr=0.01*distr/sum(distr);

tpcf=handles.Pake_t;

distr0=get_std_distr(ny2re/stretch,spc5,r);
distr0=0.01*distr0/sum(distr0);
deer=pcf2deer(distr0,kernel,r,tpcf);
% bsl=sum(deer(925:1024))/100;
% deer=deer-bsl*ones(size(deer));
deer=deer/max(deer);
sim=interp1(tpcf*dt/0.008,deer,tdip);