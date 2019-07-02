function [best_k,best_depth]=get_manual_bckg(handles)

resolution=21;

base=handles.APT_kernel;
crosstalk=handles.APT_crosstalk;
trcnorm=handles.APT_norm;
ny=handles.APT_ny;
t=handles.APT_t;

[texp,vexp,zt,phi,imo,dt]=pre_process(handles.t_orig,handles.v_orig,handles.phase,handles.imo,handles.zerotime,handles.dt);
% Determine fit range for background

t_bckg=handles.bckg_start;
t_cutoff=handles.cutoff-handles.zerotime;
ttemp=texp-t_bckg*ones(size(texp));
[tt,nofitp0]=min(abs(ttemp));
ttemp=texp-t_cutoff*ones(size(texp));
[tt,pcutoff]=min(abs(ttemp));

imagflag=get(handles.imaginary,'Value');
dcmplx=imagflag*handles.cmplx; % Complex data?

% reference deconvolution for variable-time DEER, normalization for
% constant-time DEER
if handles.ctvt,
    z=handles.A_vexp;
    ref=vexp(1,:);
    sc=max(real(ref));
    sig=vexp(2,:);
    vexp=real(sig)./real(ref)+i*imag(ref)/sc;
else
	vexp=vexp/max(real(vexp));
end;

[tmi,ztpoi]=min(abs(texp));
texp=texp(ztpoi:pcutoff);
vexp=real(vexp((ztpoi:pcutoff)));
nexp=length(texp);
texp=texp/1000;

% Adaptive background correction with depth variation

dim=handles.hom_dim;
v0=[handles.man_k handles.man_depth];
v=fminsearch(@rms_APT_bckg,v0,[],handles,texp,vexp,dim);
% v=fminsearch(@rms_Gaussian_bckg,v0,[],handles,texp,vexp,dim);
best_k=v(1);
best_depth=v(2);


% min_depth=handles.man_depth/1.5;
% max_depth=1.5*handles.man_depth;
% if max_depth>0.99,
%     max_depth=0.99;
% end;
% depths=linspace(min_depth,max_depth,resolution);
% min_k=handles.man_k/1.5;
% max_k=1.5*handles.man_k;
% ks=linspace(min_k,max_k,resolution);
% merit=zeros(resolution,resolution);
% best_k=handles.man_k;
% best_depth=handles.man_depth;
% best_merit=1.0e6;


%     for kd=1:resolution,
%         for kk=1:resolution,
%             targ=texp.^(handles.hom_dim/3);
%             bckg=(1-depths(kd))*exp(-ks(kk)*targ);
%     %         figure(13); clf;
%     %         plot(texp,vexp,'k');
%     %         hold on;
%     %         plot(texp,bckg,'r');
%             dipevo=real(vexp)-bckg;
%             dipevo=dipevo./bckg; % divide by background, eqn [13]
%             cluster=real(vexp)./bckg;
%             cluster=cluster/max(cluster);
%             [r,distr,sim]=APT(handles,texp,dipevo);
%             modsim=ones(size(sim))-sim;
%             modexp=ones(size(cluster))-cluster;
%             sc=sum(modexp.*modexp)/sum(modsim.*modexp);
%             sim=ones(size(modsim))-sc*modsim;
%             difference=sim-cluster;
%     %         figure(14); clf;
%     %         plot(texp,cluster,'k');
%     %         hold on;
%             plot(texp,sim,'r');
%             rms=sqrt(sum(difference.*difference)/(length(difference)-1));
%             merit(kd,kk)=rms;
%             if rms<best_merit,
%                 best_merit=rms;
%                 best_k=ks(kk);
%                 best_depth=depths(kd);
%             end;
%             pstr=sprintf('%s%6.3f%s%6.3f%s%6.4f','Optimizing bckg. fit parameters: k= ',ks(kk),', depth= ',depths(kd),' , r.m.s.d.: ',rms);
%             set(handles.status_line,'String',pstr);
%             % disp(pstr);
%             drawnow;
%             % keyboard
%         end;
%     end;
