%
% REGFUNCTIONAL Regularization functional constructor
%
%   fhandle = REGFUNCTIONAL('type',S,L,K,alpha)
%   Returns the cost model functional of a regularization model. The function
%   requires the signal (S), dipolar kernel (K), regularization matrix (L)
%   and regularization parameter (alpha). The type of regularization
%   functional is determined by the 'type' string argument.
%
%   fhandle = REGFUNCTIONAL('type',S,L,K,alpha,eta)
%   The Huber parameter can be specified by passing as the (eta) argument.
%
%   fhandle = REGFUNCTIONAL('type',{S1,S2,...},L,{K1,K2,...},alpha)
%   Passing multiple signals/kernels constructs minimization functional
%   as required for global fit of the regularization functionals. The global fit
%   weights are automatically computed according to their contribution
%   to ill-posedness. The calculated weights (w) can be requested as an
%   additional output argument.
%
%   fhandle = REGFUNCTIONAL('type',S,L,K,alpha,eta,w)
%   The global fit weights (w) can be manually passed to avoid computing them
%   automatically.
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function functionHandle = regfunctional(Method,Signal,RegMatrix,Kernel,RegularizationParameter,HuberParameter,weights)

if nargin<7 
   weights = []; 
end
if ~iscell(Signal)
    Signal = {Signal};
end
if ~iscell(Kernel)
    Kernel = {Kernel};
end
if nargin<6 || isempty(HuberParameter)
   HuberParameter = 1.35; 
end

%Get weights of different signals for global fitting
if isempty(weights)
weights = globalweights(Signal);
end

switch lower(Method)
    case 'tv'
        functionHandle = @(Distribution)TVfunctional(Signal,Distribution,RegMatrix,Kernel,RegularizationParameter,weights);
    case 'tikhonov'
        functionHandle = @(Distribution)TikhonovFunctional(Signal,Distribution,RegMatrix,Kernel,RegularizationParameter,weights);
    case 'huber'
        functionHandle = @(Distribution)HuberFunctional(Signal,Distribution,RegMatrix,Kernel,RegularizationParameter,HuberParameter,weights);
    otherwise
        functionHandle = Method;
end

end

%Tikhonov functional
%-----------------------------------------------------------
function [Functional,Gradient] = TikhonovFunctional(Signal,Distribution,L,Kernel,RegularizationParameter,weights)

Residual = 0;
ResidualGradient = 0;
for i = 1:length(Signal)
    Residual = Residual + weights(i)*1/2*norm(Kernel{i}*Distribution - Signal{i})^2;
    ResidualGradient =  ResidualGradient + weights(i)*Kernel{i}.'*(Kernel{i}*Distribution - Signal{i});
end

Functional = Residual + RegularizationParameter^2*1/2*norm(L*Distribution)^2;
if nargout>1
    Gradient = ResidualGradient + RegularizationParameter^2*(L')*L*Distribution;
end

end

%Total variation functional
%-----------------------------------------------------------
function [Functional,Gradient] = TVfunctional(Signal,Distribution,L,Kernel,RegularizationParameter,weights)

Residual = 0;
ResidualGradient = 0;
for i=1:length(Signal)
    Residual = Residual + weights(i)*1/2*norm(Kernel{i}*Distribution - Signal{i})^2;
    ResidualGradient = ResidualGradient + weights(i)*Kernel{i}.'*(Kernel{i}*Distribution - Signal{i});
end
Functional = Residual + RegularizationParameter^2*sum(sqrt((L*Distribution).^2 + 1e-24));
if nargout>1
    TVGradient = (L)'*diag(1./(sqrt((L*Distribution).^2 + 1e-24)))*(L)*Distribution;
    Gradient = ResidualGradient + RegularizationParameter^2*TVGradient;
end

end

%Pseudo-Huber functional
%-----------------------------------------------------------
function [Functional,Gradient] = HuberFunctional(Signal,Distribution,L,Kernel,RegularizationParameter,HuberParameter,weights)

Residual = 0;
ResidualGradient = 0;
for i=1:length(Signal)
    Residual = Residual + weights(i)*1/2*norm(Kernel{i}*Distribution - Signal{i})^2;
    ResidualGradient = ResidualGradient + weights(i)*Kernel{i}.'*(Kernel{i}*Distribution - Signal{i});
end

HuberPenalty = (sqrt(1+(L*Distribution/HuberParameter).^2)-1);
HuberFunctional = sum(HuberPenalty);
%Get the Huber norm gradient
Functional = Residual + RegularizationParameter^2*HuberFunctional;
if nargout>1
    HuberGradient =  RegularizationParameter^2*(L/(HuberParameter))'*diag(1./sqrt((L*Distribution/HuberParameter).^2 + 1))*(L/(HuberParameter))*Distribution;
    Gradient = ResidualGradient + RegularizationParameter^2*HuberGradient;
end

end