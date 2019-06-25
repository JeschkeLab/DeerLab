function deer=pcf2deer(pcf,kernel,r,t)
%
% function deer=pcf2deer(pcf,kernel,r,t)
%
% Converts the pair correlation function
% to DEER time-domain data (dipolar evolution data)
%
% see: G. Jeschke, A. Koch, U. Jonas, A. Godt,
%      J. Magn. Reson. 155, 72-82 (2002); equations [6],[10]
%
% pcf      pair correlation function G(r), the ordinate is vector r 
%          the sum over all elements of G(r)
%          is the number of spins in a spherical shell
%          with radius max(r), 
%          all elements must be smaller than 0.25
% kernel   kernel:
%          a) load pake_base40; kernel=base-ones(size(base));
%             this also provides t and r
%          b) define t and r yourself, compute kernel with make_sfm_kernel
% t        time axis of kernel, will also be time axis of DEER data
% r        distance axis of kernel, is also distance axis of
%          pair correlation function
%
% G. Jeschke, 2002
%

rsig=pcf*kernel; % multiply G(r) vector with kernel matrix
deer=exp(rsig); % eqn [10]