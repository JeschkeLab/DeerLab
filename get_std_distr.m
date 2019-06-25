function distr=get_std_distr(r0,p0,r),
%
% Computes distance distribution in standard format for pcf2deer 
% from the distribution p0(r0) that may not have equidistant sampling
% points and may not completely overlap with the range of vector r
%
% (c) G. Jeschke, 2003
%

% Standard resolution 2000 points


dr=r(2)-r(1);
rmit=min(r0);
if rmit<min(r), rmit=min(r); end;
rmat=max(r0);
if rmat>max(r), rmat=max(r); end;

rmip=round(rmit/dr);
rmap=round(rmat/dr);
rax=linspace(rmip*dr,rmap*dr,rmap-rmip+1);
distr0=interp1(r0,p0,rax,'pchip',0);
distr=0*r;
rbas=round(min(r)/dr);
distr(rmip-rbas+1:rmap-rbas+1)=distr0;



