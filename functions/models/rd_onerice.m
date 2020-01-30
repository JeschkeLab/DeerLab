function output = rd_onerice(r,param)
%
% RD_ONERICE Rician distribution parametric model
%
%   info = RD_ONERICE
%   Returns an (info) structure containing the specifics of the model.
%
%   P = RD_ONERICE(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to
%   the paramteres array (param). The required parameters can also be found
%   in the (info) structure.
%
% PARAMETERS
% name     symbol default lower bound upper bound
% --------------------------------------------------------------------------
% param(1)  <nu>    3.5     1.0         10          mean distance
% param(2)  sigma   0.7     0.1          5          standard deviation
% --------------------------------------------------------------------------
%

% This file is a part of DeerAnalysis. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.


nParam = 2;

if nargin~=0 && nargin~=2
    error('Model requires two input arguments.')
end

if nargin==0
    %If no inputs given, return info about the parametric model
    info.model  = 'Single Rice/Rician distribution';
    info.nparam  = nParam;
    info.parameters(1).name = ['Mean distance ',char(957)];
    info.parameters(1).range = [1 10];
    info.parameters(1).default = 3.5;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = ['Standard deviation ',char(963)];
    info.parameters(2).range = [0.1 5];
    info.parameters(2).default = 0.7;
    info.parameters(2).units = 'nm';
    
    output = info;
    return
end

%If user passes them, check that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

nu=param(1);
sqscale=param(2).^2;
%Compute rician/rice distribution using the zeroth order modified Bessel function of
%the first kind
P = (r./sqscale).*exp(-1/2*(r.^2 + nu.^2)./sqscale).*besseli(0,r.*nu./sqscale,1);
%The Rice distribution is zero for negative values.
P(P<0)=0;

if ~iscolumn(P)
    P = P';
end
if ~all(P==0)
    P = P/sum(P)/mean(diff(r));
end
output = P;



return