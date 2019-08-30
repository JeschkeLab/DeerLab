function output = rd_onerice(r,param)
%
% ONERICE Rician distribution parametric model
%
%   info = ONERICE
%   Returns an (info) structure containing the specifics of the model.
%
%   P = ONERICE(r,param)
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
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

nParam = 2;

if nargin==0
    %If no inputs given, return info about the parametric model
    info.Model  = 'Single Rice/Rician distribution';
    info.Equation  = ['r/',char(963),'²*exp((r² + ',char(957),'²)/(2',char(963),'²))*Bessel(r*',char(957),'/',char(963),'²)'];
    info.nParam  = nParam;
    info.parameters(1).name = ['Mean distance ',char(957)];
    info.parameters(1).range = [1 10];
    info.parameters(1).default = 3.5;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = ['Standard deviation ',char(963)];
    info.parameters(2).range = [0.1 5];
    info.parameters(2).default = 0.7;
    info.parameters(2).units = 'nm';
    
    output = info;
    
elseif nargin == 2
    
    %If user passes them, check that the number of parameters matches the model
    if length(param)~=nParam
        error('The number of input parameters does not match the number of model parameters.')
    end    
    
    nu=param(1);
    sqscale=param(2).^2;
    %Compute rician/rice distribution using the zeroth order modified Bessel function of
    %the first kind
    Distribution = (r./sqscale).*exp(-1/2*(r.^2 + nu.^2)./sqscale).*besseli(0,r.*nu./sqscale);
    %The Rice distribution is zero for negative values.
    Distribution(Distribution<0)=0;
    
    if ~iscolumn(Distribution)
        Distribution = Distribution';
    end
    if ~all(Distribution==0)
    Distribution = Distribution/sum(Distribution)/mean(diff(r));
    end
    output = Distribution;
else
    
    %Else, the user has given wrong number of inputs
    error('Model requires two input arguments.')
end


return