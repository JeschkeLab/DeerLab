%
% DD_SHELLSHELL Uniform spherical shell inside another spherical shell
%
%   info = DD_SHELLSHELL
%   Returns an (info) structure containing the specifics of the model.
%
%   P = DD_SHELLSHELL(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the paramteres array (param). The required parameters can also be found 
%   in the (info) structure.
%
% PARAMETERS
% name    symbol default lower bound upper bound
% --------------------------------------------------------------------------
% param(1)   R    1.5       0.1         20         inner shell radius
% param(2)   w1   0.5       0.1         20         1st shell thickness
% param(3)   w2   0.5       0.1         20         2nd shell thickness
% --------------------------------------------------------------------------
%
%   See: D.R. Kattnig, D. Hinderberger, Journal of Magnetic Resonance, 230 (2013), 50-63 
%        http://doi.org/10.1016/j.jmr.2013.01.007
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function output = dd_shellshell(r,param)

nParam = 3;

if nargin~=0 && nargin~=2 
    error('Model requires two input arguments.')
end

if nargin==0
    %If no inputs given, return info about the parametric model
    info.model  = 'Uniform spherical shell inside another spherical shell';
    info.nparam  = nParam;
    info.parameters(1).name = 'Sphere radius R';
    info.parameters(1).range = [0.1 20];
    info.parameters(1).default = 1.5;
    info.parameters(1).units = 'nm';

    info.parameters(2).name = '1st Shell thickness';
    info.parameters(2).range = [0.1 20];
    info.parameters(2).default = 0.5;
    info.parameters(2).units = 'nm';
    
    info.parameters(3).name = '2nd Shell thickness';
    info.parameters(4).range = [0.5 20];
    info.parameters(4).default = 0.5;
    info.parameters(4).units = 'nm';
    output = info;
    return
end

% Assert that the number of parameters matches the model
if length(param)~=nParam
  error('The number of input parameters does not match the number of model parameters.')
end

%Parse input
validateattributes(r,{'numeric'},{'nonnegative','increasing','nonempty'},mfilename,'r')

% Compute the model distance distribution
R1 = param(1);
w1 = param(2);
w2 = param(3);

R2 = R1 + w1;
R3 = R2 + w2;

delta21 = R2^3 - R1^3;
q21 = delta21*pbs(r,R1,R2);
delta31 = R3^3 - R1^3;
q31 = delta31*pbs(r,R1,R3);
delta32 = R3^3 - R2^3;
q32 = delta32*pbs(r,R2,R3);

P = R1^3*q21 - R1^3*q31 + R2^3*q32;
P = P/(delta21*delta32);

if ~all(P==0)
P = P/trapz(r,P);
end

output = P;

function P = pbs(r,R1,R2)

P = zeros(numel(r),1);
%Case1
idx = r >= 0 & r < min(2*R1,R2 - R1); 
P(idx) = 12*r(idx).^3*R1^2 - r(idx).^5;

%Case2
idx = r >= R2 - R1 & r < 2*R1;
P(idx) = 8*r(idx).^2*(R2^3 - R1^3) - 3*r(idx)*(R2^2 - R1^2)^2 - 6*r(idx).^3*(R2 - R1)*(R2 + R1);

%Case3
idx = r >= 2*R1 & r < R2 - R1;
P(idx) = 16*r(idx).^2*R1^3;

%Case4
idx = r >= max(R2 - R1,2*R1) & r < R1 + R2;
P(idx) = r(idx).^5 - 6*r(idx).^3*(R2^2 + R1^2) + 8*r(idx).^2*(R2^3 + R1^3) - 3*r(idx)*(R2^2 - R1^2)^2;

P = P*3/(16*R1^3*(R2^3 - R1^3));

end

end
