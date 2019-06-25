function rms=rms_APT_bckg(v,handles,texp,vexp,hom_dim)
% computes the r.m.s.d. of the form factor computed from an APT-derived
% distance distribution from the experimental form factor
% a homogeneous background correction with dimension dim is applied, whose
% parameters are given in vector v
%
% v(1)  decay time constant
% v(2)  modulation depth
% texp  time axis of zero-time corrected, normalized and cutoff original data
% vexp  zero-time corrected and cutoff original data
%
% (c) G. Jeschke, 2008

if v(1)<0.001,
    rms=1e6;
    return;
end;

if v(2)<0.01 || v(2)>0.99,
    rms=1e6;
    return;
end;

v1(1)=v(1);
v1(2)=1;
bckg=decaynD(v1,texp,hom_dim);
bckg=(1-v(2))*bckg;

% targ=texp.^(dim/3);
% bckg=(1-v(2))*exp(-v(1)*targ);
%         figure(13); clf;
%         plot(texp,vexp,'k');
%         hold on;
%         plot(texp,bckg,'r');
dipevo=real(vexp)-bckg;
dipevo=dipevo./bckg; % divide by background, eqn [13]
cluster=real(vexp)./bckg;
cluster=cluster/max(cluster);
[r,distr,sim]=APT(handles,texp,dipevo);
modsim=ones(size(sim))-sim;
modexp=ones(size(cluster))-cluster;
sc=sum(modexp.*modexp)/sum(modsim.*modexp);
sim=ones(size(modsim))-sc*modsim;
difference=sim-cluster;
%         figure(14); clf;
%         plot(texp,cluster,'k');
%         hold on;
%         plot(texp,sim,'r');
rms=sqrt(sum(difference.*difference)/(length(difference)-1));
pstr=sprintf('%s%6.3f%s%6.3f%s%6.4f','Optimizing bckg. fit parameters: k= ',v(1),', depth= ',v(2),' , r.m.s.d.: ',rms);
set(handles.status_line,'String',pstr);
% disp(pstr);
drawnow;
% keyboard
