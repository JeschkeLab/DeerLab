function [texp,vexp,zt,phi,imo,dt,renorm]=pre_process(t0,v0,phi,imo,zt,dt)
%
% Preprocessing of DEER data
% 
% pre-processed data are phase-corrected, shifted on the time axis to make
% the maximum coincide with time zero, have a time increment that is a
% positive integer multiple of 4 ns, and have a maximum of 2048 data points
%
% when data reduction is required, information from all primary data points
% is included by smoothing, if requested data points fall in between
% primary data points, spline interpolation is used, data are normalized to
% the maximum of the real part after phase correction
%
% Input parameters:
% t0    time axis of original experimental data
% v0    original experimental data (echo integral)
% phi   phase shift (optional), if not given, an automatic phase correction
%       is performed that minimizes the imaginary part
% imo   (optional) offset to be subtracted from imaginary part before phase 
%       correction, if not given, it is determined automatically
% zt    zero time shift (optional), if not given, an automatic
%       determination is performed that minimizes the first moment around
%       the maximum
% dt    time increment of the pre-processed data, must be a positive
%       multiple of eight (optional), defaults to the lowest positive
%       multiple of eight that leads to a data set with at most 2048 points
% renorm    renormalization factor for maximum determination with
%           polynomial smoothing
%
% Additional output parameters:
% texp  zero-shift corrected experimental time axis
% vexp  phase corrected data

max_points = 2048; % maximum number of data points
t_grain = 4; % minimum time increment, should coincide with handles.time_grain in DeerAnalysis.m

add_noise = 0; % tweak for adding noise on loading, only for test purposes

fit_offset=false;

[m,n]=size(v0); % to distinguish between constant- and variable-time DEER

if m>1 
    vx=v0(1,:); 
else
    vx=v0;
end

% get phase value, if not provided
if nargin<3
    % phi=fminbnd(@rms_phi,-pi/2+1e-4,pi/2,[],vx);
    last=v0(length(v0));
    phi0=atan2(imag(last),real(last));
    im_off_0=0;
    v=[phi0,im_off_0];
    pstart=round(length(vx)/8); % use only last 7/8 of data for phase/offset correction
    if ~fit_offset
        v=v(1);
    end
    v=fminsearch(@rmsd_phi_offset,v,[],vx(pstart:end));
    phi=v(1);
    if fit_offset
        imo=v(2);
    else
        imo=0;
    end
    if sum(real(vx*exp(1i*phi)))<0, phi=phi+pi; end
end
% make phase correction and normalize
if m>1
    ref=(v0(1,:)-1i*imo)*exp(1i*phi);
    sc=max(real(ref));
    sig=(v0(2,:)-1i*imo)*exp(1i*phi);
    vexp=real(sig)./real(ref)+1i*imag(ref)/sc;
else
	vexp=(v0-1i*imo)*exp(1i*phi);
    sc = 1/max(real(vexp));
%     phi=fminbnd(@rms_phi,-pi/2+1e-4,pi/2,[],vx); % old version!
% 	vexp=v0*exp(i*phi);
%     imo=0;
	vexp=sc*vexp;
end

t1=min(t0):1:max(t0);
v1=interp1(t0,real(vexp),t1,'spline',real(vexp(1)));
% get zero time, if not provided
expdat=real(v1);
if nargin<5
    % Determine maximum
    [~,mp]=max(expdat);
    nmp=1;
    % Determine time zero by moment analysis
    if mp>1 && mp<length(expdat)
        dmi=mp-1;
        dma=length(expdat)-mp;
        dm=dmi; if dma<dm, dm=dma; end
        maxd=floor(dm/2);
        dev=1e20;
        nmp=mp;
        for k=-maxd+mp:maxd+mp
            summ=0;
            for l=-maxd:maxd
                summ=summ+expdat(k+l)*l;
            end
            if abs(summ)<dev, dev=abs(summ); nmp=k; end
        end
    end
    zt=t1(nmp);
else
    [~,nmp] = min(abs(t1-zt));
end

renorm = 1;
% get smoothed estimate of maximum
if nmp > 10
    xax = 1:2*nmp;
    [p,~] = polyfit(xax,expdat(1:2*nmp),11);
    smoothed = polyval(p,xax);
    renorm = smoothed(nmp);
    vexp = vexp/smoothed(nmp);
%     expdat = expdat/smoothed(nmp);
%     smoothed = smoothed/smoothed(nmp);
%     figure(13); clf;
%     plot(expdat,'k');
%     hold on;
%     plot([nmp,nmp],[0,1],'b');
%     plot(1:2*nmp,smoothed,'r');
end



% Correct time axis
t0=t0-zt*ones(size(t0));
ta=t0(1); % real start time
te=t0(length(t0)); % real end time

dt0=t0(2)-t0(1);

% get dt if not provided
if nargin<6
    dt=dt0;
    if dt<t_grain, dt=t_grain; end
    if mod(dt,t_grain)~=0
        dt=t_grain*ceil(dt/t_grain);
    end
end
redfac0=dt/dt0; % data reduction factor for requested dt

max_nexp=max_points; % no more than the specified maximum number of data points
if n>max_nexp
    redfac=ceil(n/max_nexp);
    if redfac>redfac0 % requested dt is too small
        dt=redfac*dt0;
    end
    if mod(dt,t_grain)~=0 % next multiple of eight
        dt=t_grain*ceil(dt/t_grain);
    end
    redfac=round(dt/dt0);
    % Data reduction with smoothing
    pas=pascal(redfac);
    b=zeros(1,redfac);
    for k=1:redfac
        b(k)=pas(redfac-k+1,k);
    end
    b=b/sum(b);
    z0=vexp;
    a=1;
    [m2,n2]=size(z0);
    z=zeros(m2,n2);
    for kk=1:m2
        z0a=z0(kk,:);
        z0a=filter(b,a,z0a);
        z(kk,:)=z0a;
        z(kk,1)=z0(kk,1);
    end
else
    z=vexp;
end

% Make new time axis
mint=dt*ceil(ta/dt);
if mint<ta, mint=ta; end
maxt=dt*floor(te/dt);
% while maxt>te, maxt=maxt-dt; end;
if maxt>te, maxt=te; end
newn=1+floor((maxt-mint)/dt);
maxt=mint+dt*(newn-1);
texp=mint:dt:maxt;

% Interpolation onto new time axis
vexp=zeros(m,length(texp));
[m2,~]=size(z);
for kk=1:m2
    vexp(kk,:)=interp1(t0,z(kk,:),texp,'spline',z(kk,1));
end

vexp = vexp + add_noise*randn(size(vexp));

function rmsd=rmsd_phi_offset(v,tr)
% Computes root mean square deviation of the imaginary part of
% phase-corrected data from zero
% before phase correction, an offset is subtracted from the imaginary part
%
% v(1)  phase correction phi (rad)
% v(2)  offset
% tr    complex data trace
%
% rmsd  root mean square deviation of imaginary part from zero
%
% G. Jeschke, 2009

if length(v)>1
    tr=tr-1i*v(2);
    itr=imag(tr*exp(1i*v(1)));
    rmsd=sqrt(sum(itr.*itr));
else
    itr=imag(tr*exp(1i*v(1)));
    rmsd=sqrt(sum(itr.*itr));
end
