function sim=get_td_fit(handles,r0,distr0),

if length(r0)*length(distr0)==0,
    sim=[];
    return;
end;

r=handles.Pake_r;
kernel=handles.Pake_kernel;
tpcf=handles.Pake_t;
tdip=handles.A_tdip;
dt=tdip(2)-tdip(1); % time increment
dt0=tpcf(2)-tpcf(1);
stretch=(dt/dt0)^(1/3);

% disp(size(r0));
% disp(size(stretch));
% disp(size(distr0));
% disp(size(r));
% disp('MIau');
distr=get_std_distr(r0/stretch,distr0,r);
distr=0.01*distr/sum(distr);


deer=pcf2deer(distr,kernel,r,tpcf);
% bsl=sum(deer(925:1024))/100;
% deer=deer-bsl*ones(size(deer));
% deer=deer/max(deer);
sim=interp1(tpcf*dt/0.008,deer,tdip);

