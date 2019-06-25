function rms=rms_Gauss_bckg(v,handles,texp,vexp,dim)
% computes the r.m.s.d. of the form factor computed from a Tikhonov-derived
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

if v(2)<0.01 | v(2)>0.99,
    rms=1e6;
    return;
end;

targ=texp.^(dim/3);
bckg=(1-v(2))*exp(-v(1)*targ);
%         figure(13); clf;
%         plot(texp,vexp,'k');
%         hold on;
%         plot(texp,bckg,'r');
dipevo=real(vexp)-bckg;
dipevo=dipevo./bckg; % divide by background, eqn [13]
cluster=real(vexp)./bckg;
cluster=cluster/max(cluster);
v0=[handles.rmin,handles.sr];
[v1,rms]=fminsearch(@rms_Gaussian,v0,[],handles,cluster);
pstr=sprintf('%s%6.3f%s%6.3f%s%6.4f','Optimizing bckg. fit parameters: k= ',v(1),', depth= ',v(2),' , r.m.s.d.: ',rms);
set(handles.status_line,'String',pstr);
drawnow;

function rms=rms_Gaussian(v,handles,cluster)
%
r=linspace(1.5,8,201);
gauss0=(r-v(1)*ones(size(r)))/v(2);
distr=exp(-gauss0.^2);
distr=0.01*distr/sum(distr); % normalization

% Simulate fitted dipolar evolution function
sim=get_td_fit(handles,r,distr);

modsim=ones(size(sim))-sim;
modexp=ones(size(cluster))-cluster;
sc=sum(modexp.*modexp)/sum(modsim.*modexp);
if sc<0, rms=1e6; return; end;
sim=ones(size(modsim))-sc*modsim;
diff=sim-cluster;
rms=sqrt(sum(diff.*diff))/(length(diff)-1);
