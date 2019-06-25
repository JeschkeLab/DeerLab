function filtered=long_pass_filter(handles,tdip,dipevo);
% Low-pass filtering to suppress proton modulations, suppresses distances
% below threshold r_filter_min (default: 1.5 nm)
%
% tdip    time axis
% dipevo  dipolar evolution function
%

rmin=handles.longpass_min;
dipevo0=dipevo;

if mod(length(dipevo),2)==1,
    dipevo=dipevo(1:length(dipevo)-1);
    tdip=tdip(1:length(tdip)-1);
end;

cutoff=52.04/(rmin^3); % compute cutoff frequency in MHz (low pass) from minimum distance

dt=tdip(2)-tdip(1);
ny=linspace(-1/(2*dt),1/(2*dt),length(tdip));
	
ny1=ny+cutoff*ones(size(ny));
[voi,filta]=min(abs(ny1));

ny2=ny-cutoff*ones(size(ny));
[voi,filte]=min(abs(ny2));
		
	
spc=fftshift(fft(dipevo));
	
spc1=spc;
frange=spc1(1:filta-1);
fx=linspace(1,length(frange),length(frange));
[p,s]=polyfit(fx,frange,3);
corr=polyval(p,fx);
spc1(1:filta-1)=zeros(size(spc1(1:filta-1)));
spc1(1:filta-1)=corr;
frange=spc1(filte:length(spc1));
fx=linspace(1,length(frange),length(frange));
[p,s]=polyfit(fx,frange,3);
corr=polyval(p,fx);
spc1(filte:length(spc1))=zeros(size(spc1(filte:length(spc1))));
spc1(filte:length(spc1))=corr;

% figure(13); clf;
% plot(real(spc),'k');
% hold on;
% plot(imag(spc),'b');
% plot(real(spc1),'r');
% plot(imag(spc1),'m');
%keyboard;

fid=ifft(fftshift(spc1));
filtered=real(fid);	

if length(dipevo0)>length(filtered),
    filtered=[filtered dipevo0(length(dipevo0))];
end;
