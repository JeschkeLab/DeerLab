function output = rd_tworice(r,param)
%
% TWORICE Sum of two rician distributions parametric model
%
%   info = TWORICE
%   Returns an (info) structure containing the specifics of the model.
%
%   P = TWORICE(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to
%   the paramteres array (param). The required parameters can also be found
%   in the (info) structure.
%
% PARAMETERS
% name      symbol default lower bound upper bound
% --------------------------------------------------------------------------
% param(1)  <nu1>    2.5     1.0        10         mean distance
% param(2)  sigma1   0.4     0.1        5          standard deviation
% param(3)  <nu2>    4.0     1.0        10         mean distance
% param(4)  sigma2   0.4     0.1        5          standard deviation
% param(5)  p1       0.5     0          1          fraction of pairs at 1st distance
% --------------------------------------------------------------------------
%

% This file is a part of DeerAnalysis. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.


nParam = 5;

if nargin~=0 && nargin~=2
    error('Model requires two input arguments.')
end

if nargin==0
    %If no inputs given, return info about the parametric model
    info.model  = 'Two Rice/Rician distributions';
    info.Equation  = ['A*r/',char(963),'1²*exp((r²+',char(957),'1²)/(2',char(963),'1²))*Bessel(r*',char(957),'1/',char(963),'1²)'...
        '+ (1-A)*r/',char(963),'2²*exp((r² + ',char(957),'2²)/(2',char(963),'2²))*Bessel(r*',char(957),'2/',char(963),'2²)' ];
    info.nparam  = nParam;
    info.parameters(1).name = ['Mean distance ',char(957),'1 1st Rician'];
    info.parameters(1).range = [1 10];
    info.parameters(1).default = 2.5;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = ['Standard deviation ',char(963),'1 1st Rician'];
    info.parameters(2).range = [0.1 5];
    info.parameters(2).default = 0.7;
    info.parameters(2).units = 'nm';
    
    info.parameters(3).name = ['Mean distance ',char(957),'2 2nd Rician'];
    info.parameters(3).range = [1 10];
    info.parameters(3).default = 4.0;
    info.parameters(3).units = 'nm';
    
    info.parameters(4).name = ['Standard deviation ',char(963),'2 2nd Rician'];
    info.parameters(4).range = [0.1 5];
    info.parameters(4).default = 0.7;
    info.parameters(4).units = 'nm';
    
    info.parameters(5).name = 'Relative amplitude A 1st Rician';
    info.parameters(5).range = [0 1];
    info.parameters(5).default = 0.5;
    
    output = info;
    return;
end

%If user passes them, check that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

nu=param(1);
sqscale=param(2).^2;
%Compute rician/rice distribution using the zeroth order modified Bessel function of
%the first kind
Rician1 = (r./sqscale).*exp(-1/2*(r.^2 + nu.^2)./sqscale).*besseli(0,r.*nu./sqscale);
%The Rice distribution is zero for negative values.
Rician1(Rician1<0)=0;

nu=param(3);
sqscale=param(4).^2;
%Compute rician/rice distribution using the zeroth order modified Bessel function of
%the first kind
Rician2 = (r./sqscale).*exp(-1/2*(r.^2 + nu.^2)./sqscale).*besseli(0,r.*nu./sqscale);
%The Rice distribution is zero for negative values.
Rician2(Rician2<0)=0;

%Construct distance distribution
Distribution = param(5)*Rician1 + max(1-param(5),0)*Rician2;

if ~iscolumn(Distribution)
    Distribution = Distribution';
end
if ~all(Distribution==0)
    Distribution = Distribution/sum(Distribution)/mean(diff(r));
end
output = Distribution;

return