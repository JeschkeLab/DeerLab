%
% RD_SPHEREPOINT One particle distanced from particles distributed on a
% sphere 
%
%   info = RD_SPHEREPOINT
%   Returns an (info) structure containing the specifics of the model.
%
%   P = RD_SPHEREPOINT(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the paramteres array (param). The required parameters can also be found 
%   in the (info) structure.
%
% PARAMETERS
% name    symbol default lower bound upper bound
% --------------------------------------------------------------------------
% param(1)   R    1.5       0.1         20         sphere radius
% param(2)   d    3.5       0.1         20         distance from sphere center to point
% --------------------------------------------------------------------------
%
%   See: D.R. Kattnig, D. Hinderberger, Journal of Magnetic Resonance, 230 (2013), 50-63 
%        http://doi.org/10.1016/j.jmr.2013.01.007
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.


function output = rd_spherepoint(r,param)

nParam = 2;

if nargin~=0 && nargin~=2 
    error('Model requires two input arguments.')
end

if nargin==0
    %If no inputs given, return info about the parametric model
    info.model  = 'Distance between two points inside/outside a sphere';
    info.nparam  = nParam;
    info.parameters(1).name = 'Sphere radius R';
    info.parameters(1).range = [0.1 20];
    info.parameters(1).default = 1.5;
    info.parameters(1).units = 'nm';

    info.parameters(2).name = 'Distance from sphere center to point';
    info.parameters(2).range = [0.1 20];
    info.parameters(2).default = 3.5;
    info.parameters(2).units = 'nm';
    output = info;
    return
end

% Assert that the number of parameters matches the model
if length(param)~=nParam
  error('The number of input parameters does not match the number of model parameters.')
end

% Compute the model distance distribution
R = param(1);
d = param(2);
P = zeros(numel(r),1);
idx = r >= d - R & r<= d + R; 
P(idx) = 3*r(idx).*(R^2 - (d - r(idx)).^2)./(4*d*R.^3);

if ~all(P==0)
P = P/sum(P)/mean(diff(r));
end

output = P;

return