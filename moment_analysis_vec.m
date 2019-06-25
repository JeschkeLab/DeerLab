function momvec=moment_analysis(x0,y0),
%
% Computes the first three moments m1, m2, and m3 of distribution y(x)
% also determines linewidths lw50 in which 50% and lw90 in which 90%
% of all intensity is contained
%
% (c) G. Jeschke, 2003
%

% Standard resolution 2000 points

x=linspace(min(x0),max(x0),2000);
y=interp1(x0,y0,x,'pchip');

lw50=0;
lw90=0;

% Compute the moments
m1=sum(y.*x)/sum(y);
delx=x-m1*ones(size(x));
m2=sum(y.*delx.*delx)/sum(y);
m3=sum(y.*delx.*delx.*delx)/sum(y);

% Approximate linewidth is square root of second moment
if m2>0,
    lw0=sqrt(m2);
else,
    lw0=0.5;
end;
dx=lw0/100; % interpolate data to intervals of 1% of approximate linewidth
np=round((max(x)-min(x))/dx);
refx=linspace(min(x),max(x),np);
refy=interp1(x,y,refx,'pchip');
refy=refy/sum(refy); % normalize distribution to unity integral
drefx=refx-m1*ones(size(refx));
[mix,midpoint]=min(abs(drefx));
symline=zeros(size(refx));
lleft=fliplr(refy(1:midpoint-1));
lright=refy(midpoint+1:length(refy));
symline(1:length(lleft))=lleft;
symline(1:length(lright))=symline(1:length(lright))+lright;
symintg=cumsum(symline);
det50=symintg-0.5*ones(size(symintg));
det90=symintg-0.9*ones(size(symintg));
[mi50,p50]=min(abs(det50));
lw50=2*dx*p50;
[mi90,p90]=min(abs(det90));

lw90=2*dx*p90;

mi=min(refy);
ma=max(refy);
yy=[mi-0.1*(ma-mi) ma+0.1*(ma-mi)];

% figure(1); clf;
% plot(refx,refy,'k');
% axis([min(refx),max(refx),mi-0.2*(ma-mi),ma+0.2*(ma-mi)]);
% axis([1,8,mi-0.2*(ma-mi),ma+0.2*(ma-mi)]);
% hold on
% plot([m1 m1],yy,'r');
% plot([m1-lw50/2 m1-lw50/2],yy,'b:');
% plot([m1+lw50/2 m1+lw50/2],yy,'b:');
% plot([m1-lw90/2 m1-lw90/2],yy,'g:');
% plot([m1+lw90/2 m1+lw90/2],yy,'g:');

momvec=[m1,m2,m3,lw50,lw90];




