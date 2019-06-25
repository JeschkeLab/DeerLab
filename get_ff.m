function sim=get_ff(r0,distr0,t,kernel,r,tpcf)

if ~exist('r','var'),
    r = r0;
end;

if ~exist('tpcf','var'),
    tpcf = t;
end;

dt=t(2)-t(1); % time increment
dt0=tpcf(2)-tpcf(1);
stretch=(dt/dt0)^(1/3);

distr=get_TK_std_distr(r0/stretch,distr0,r);
distr=0.01*distr/sum(distr);


rsig=distr*kernel; % multiply G(r) vector with kernel matrix
deer=exp(rsig); % eqn [10]

sim=interp1(tpcf*dt/dt0,deer,t);
sim = sim - 1;
sim = sim/max(sim);

