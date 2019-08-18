%
% REGOPERATOR Compute discrete derivative regularization operators 
% 
%   L = REGOPERATOR(n,d)
%   Computes the discrete approximation L to the derivative operator 
%   of order d on a regular grid with n points, i.e. L is (n-d)-by-n. 
% 
%   [L,W] = REGOPERATOR(n,d)
%   If requested, also computes W, an orthonormal basis for the null 
%   space of L. 
%
% Adapted from Christian Hansen
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function [RegMatrix,OrthNormBase] = regoperator(Dimension,Order) 

 
% Initialization. 
if Order<0
    error ('Order d must be nonnegative')
end 
 
% Zero'th derivative. 
if Order==0
    RegMatrix = speye(Dimension);
    OrthNormBase = zeros(Dimension,0); 
    return 
end 
 
% Compute L. 
c = [-1,1,zeros(1,Order-1)]; 
nd = Dimension-Order; 
for i=2:Order, c = [0,c(1:Order)] - [c(1:Order),0]; end 
RegMatrix = sparse(nd,Dimension); 
for i=1:Order+1 
  RegMatrix = RegMatrix + sparse(1:nd,[1:nd]+i-1,c(i)*ones(1,nd),nd,Dimension); 
end 
 
% If required, compute the null vectors W via modified Gram-Schmidt. 
if (nargout==2) 
  OrthNormBase = zeros(Dimension,Order); 
  OrthNormBase(:,1) = ones(Dimension,1); 
  for i=2:Order, OrthNormBase(:,i) = OrthNormBase(:,i-1).*[1:Dimension]'; end 
  for k=1:Order 
     OrthNormBase(:,k) = OrthNormBase(:,k)/norm(OrthNormBase(:,k)); 
     OrthNormBase(:,k+1:Order) = OrthNormBase(:,k+1:Order) - OrthNormBase(:,k)*(OrthNormBase(:,k)'*OrthNormBase(:,k+1:Order)); 
  end 
end 
