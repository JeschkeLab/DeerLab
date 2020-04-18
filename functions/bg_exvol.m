%
% BG_EXVOL Excluded-volume model
%
%   info = BG_EXVOL
%   Returns an (info) structure containing the specifics of the model.
%
%   B = BG_EXVOL(t,param)
%   Computes the N-point model (B) from the N-point time axis (t) according to
%   the paramteres array (param). The required parameters can also be found
%   in the (info) structure.
%
% PARAMETERS
% name    symbol  default lower bound upper bound
% ----------------------------------------------------------------------------
% PARAM(1)  R       1         0.1         20        distance of closest approach (nm)
% PARAM(1)  lamc    50        0.01       1000       lambda*concentration (uM)
% ----------------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function output = bg_exvol(t,param)

nParam = 2;

if nargin~=0 && nargin~=2
    error('Model requires two input arguments.')
end

if nargin==0
    % If no inputs given, return info about the parametric model
    info.model  = 'Excluded volume';
    info.nparam  = nParam;
    
    info.parameters(1).name = 'distance of closest approach R';
    info.parameters(1).range = [0.1 20];
    info.parameters(1).default = 1;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = 'concentration of excited spins';
    info.parameters(2).range = [0.01 1000];
    info.parameters(2).default = 50;
    info.parameters(2).units = 'uM';
    
    output = info;
    return
end

% If user passes them, check that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters (%d) does not match the number of model parameters (%d).',...
        length(param),nParam)
end

% Load precalculated reduction factor look-up table (Kattnig Eq.(18))
% To regenerate look-up table, use private/bg_exvol_alpha
persistent exvol
if isempty(exvol)
    load('bg_exvol','exvol');
end

% Get parameters
R = param(1); % nm
lambda_c = param(2); % uM

NA = 6.02214076e23; % Avogadro constant, mol^-1
lambda_c = lambda_c*1e-6*1e3*NA; % umol/L -> mol/L -> mol/m^3 -> spins/m^3

A = (mu0/4/pi)*(gfree*bmagn)^2/hbar; % Eq.(6); m^3 rad/s

% Calculate reduction factor (Eq.(18))
if R==0
    alpha = 1;
else
    dR = A*abs(t*1e-6)/(R*1e-9)^3; % unitless
    
    % Use interpolation of look-up table for small dR
    small = dR<max(exvol.dR);
    alpha = zeros(size(dR));
    alpha(small) = interp1(exvol.dR,exvol.alpha,dR(small),'makima');
    
    % For large dR, use limiting dR->inf expression
    alpha(~small) = 1 - (3/2/pi)*sqrt(3)./dR(~small);
end

K = 8*pi^2/9/sqrt(3)*A*abs(t*1e-6).*alpha; % Eq.(17)
V = exp(-lambda_c*K); % Eq.(13)

if ~iscolumn(V)
    V = V.';
end

output = V;

return
