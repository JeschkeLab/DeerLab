function [base,tnorm,ny,t,crosstalk]=getAPTkernel(numdat)
% Computes kernel for approximate Pake transformation
%
% numdat    number of data points in time domain
%
% base      kernel data
% norm      normalization constants, eqn [19] of Ref. (1)
% ny        dipolar frequency axis (MHz)
% t         time axis (µs)
% crosstalk cross-talk matrix, eqn [20] of Ref. 1
%
% (1) G. Jeschke, A. Koch, U. Jonas, A. Godt, J. Magn. Reson. 155, 72-82 (2002)
%
% (c) G. Jeschke, 2001,2019
%

t=linspace(0,(numdat-1)*0.008,numdat); % time axis, increment 8 ns
dny=1/(2*max(t)); % eqn [23] of Ref. (1)
numny=floor(numdat/2)-2; % length of dipolar frequency axis
ny=linspace(1,numny,numny);
ny=dny*(ny+1/4*ones(1,numny)); % still eqn [23]
tnorm=zeros(1,numny); % initialize vector of normalization constants
base=zeros(numny,numdat); % initialize kernel array
for k=1:numny % loop over dipolar frequency values
   wdd=2*pi*ny(k); % angular frequency
   if mod(k,10) == 0
   end
   for x=0:0.001:1 % loop over theta values
      ww=wdd*(3*x^2-1); % current dipolar frequency
      base(k,:)=base(k,:) + cos(ww*t); % add data trace to kernel
   end
end
for k=1:numny % normalize kernel traces to value at time origin
   base(k,:)=base(k,:)/base(k,1);
   tnorm(k)=sum(base(k,:).*base(k,:).*t); % compute normalization constant, eqn [19]
end
[m,~]=size(base); % size of kernel
crosstalk=zeros(m,m); % initialize crosstalk matrix
for k=1:m % compute crosstalk matrix, eqn [20]
    for l=1:m
        mu=base(k,:);
        crosstalk(k,l)= sum(mu.*base(l,:).*t)/tnorm(k);
    end
end
condition=cond(crosstalk); % condition number of crosstalk matrix, currently unused
