function deer=get_DEER(r0,distr,moddepth,kernel,r,t)
%
% Computes a DEER form factor with given modulation depth moddepth from a
% distance distribution distr with corresponding distance axis r0, the
% kernel table kernel and the corresponding distance and time axes r and t
% must be provided
%
% (c) G. Jeschke, 2008

distr_1=get_std_distr(r0,distr,r);
distr_1=0.01*distr_1/sum(distr_1);
deer0=pcf2deer(distr_1,kernel,r,t);

moddeer=ones(size(deer0))-deer0;
deer=ones(size(moddeer))-moddepth*moddeer/0.01;
