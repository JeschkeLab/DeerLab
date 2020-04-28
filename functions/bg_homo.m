%
% BG_HOMO Multi-pulse DEER background in a homogenous medium
%
%   info = BG_HOMO
%   Returns an (info) structure containing the specifics of the model.
%
%   B = BG_HOMO(t,param)
%   B = BG_HOMO(t,param,lambda)
%   Computes the N-point model (B) from the N-point time axis (t) according to
%   the paramteres array (param). The required parameters can also be found
%   in the (info) structure. The pathway amplitude (lambda) can be
%   included, if not given the default lambda=1 will be used.
%
% PARAMETERS
% name    symbol default lower bound upper bound
% ------------------------------------------------------------------
% PARAM(1)  c       50        0.01       5000     concentration (uM)
% ------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.



function output = bg_homo(t,param,lambda)

nParam = 1;

if all(nargin~=[0 2 3])
    error('Model requires at least two input arguments.')
end

if nargin==0
    %If no inputs given, return info about the parametric model
    info.model  = 'Multi-pulse DEER background in a homogenous medium';
    info.nparam  = nParam;
    info.parameters(1).name = 'Spin concentration';
    info.parameters(1).range = [0.01 5000];
    info.parameters(1).default = 50;
    info.parameters(1).units = 'uM';
    
    output = info;
    return
end

if nargin<3
    lambda = 1;
end

% If user passes them, check that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

% If necessary inputs given, compute the model background
conc = param(1); %uM
NA = 6.02214076e23; % Avogadro constant, mol^-1
conc = conc*1e-6*1e3*NA; % umol/L -> mol/L -> mol/m^3 -> spins/m^3

muB = 9.2740100783e-24; % Bohr magneton, J/T (CODATA 2018 value);
mu0 = 1.25663706212e-6; % magnetic constant, N A^-2 = T^2 m^3 J^-1 (CODATA 2018)
h = 6.62607015e-34; % Planck constant, J/Hz (CODATA 2018)
ge = 2.00231930436256; % free-electron g factor (CODATA 2018 value)

D = (mu0/4/pi)*(muB*ge)^2/h*1e-6; %m^3 mus^-1

% Compute constant
kappa = 8*pi^2/9/sqrt(3)*D;

% Compute background function
B = exp(-lambda*conc*kappa*abs(t));
B = B(:);
output = B;

return