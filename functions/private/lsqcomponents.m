%
% LSQCOMPONENTS Construct input arguments for NNLS solvers
%
%   [Q,KtS] = LSQCOMPONENTS(S,K,L,alpha,'type')
%   Computes the components required by non-negative least-squares (NNLS)
%   solvers (Q) and (KtS) for fitting a regularization model. The function
%   requires the signal (S), dipolar kernel (K), regularization matrix (L)
%   and regularization parameter (alpha). The type of regularization
%   functional is determined by the 'type' string argument.
%
%   [Q,KtS] = LSQCOMPONENTS(S,K,L,alpha,'type',eta)
%   The Huber parameter can be specified by passing as the (eta) argument.
%
%   [Q,KtS,w] = LSQCOMPONENTS({S1,S2,...},{K1,K2,...},L,alpha,'type',eta)
%   Passing multiple signals/kernels constructs the LSQ components (Q) and (KtS)
%   as required for global fit of the regularization functionals. The global fit
%   weights are automatically computed according to their contribution
%   to ill-posedness. The calculated weights (w) can be requested as an
%   additional output argument.
%
%   [Q,KtS] = LSQCOMPONENTS(S,K,L,alpha,'type',eta,w)
%   The global fit weights (w) can be manually passed to avoid computing them
%   automatically.
%

% This file is a part of DeerAnalysis. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.


function [Q,KtS,weights] = lsqcomponents(Signal,Kernel,RegMatrix,RegParam,RegType,HuberParam,weights)

%Ensure that signals and kernel are in a cell array
if ~iscell(Signal)
    Signal = {Signal};
end
if ~iscell(Kernel)
    Kernel = {Kernel};
end
%If Huber parameter not given, just use the default
if nargin<6 || isempty(HuberParam)
   HuberParam = 1.35; 
end
%If Huber parameter not given, just use the default
if nargin<7 
    weights = [];
end
%Prepare
nSignals = length(Signal);
distDim = length(RegMatrix);
KtS = zeros(distDim,1);
GramMatrix = zeros(distDim,distDim);

%Get weights of different signals for global fitting
if isempty(weights)
weights = globalweights(Signal);
end
%Compute the terms depending on different signals
for i=1:nSignals
    KtS = KtS + weights(i)*Kernel{i}.'*Signal{i};
    GramMatrix = GramMatrix + weights(i)*Kernel{i}.'*Kernel{i};
end

%Compute then the LSQ components needed by NNLS optimizers
switch lower(RegType)
    
    case 'tikhonov'
        Q = GramMatrix + RegParam^2*(RegMatrix.'*RegMatrix);
        
    case 'tv'
        localDistribution = zeros(distDim,1);
        for j=1:500
            prev = localDistribution;
            %Compute pseudoinverse and unconst. distribution recursively
            TVterm = RegMatrix.'*((RegMatrix./sqrt((RegMatrix*localDistribution).^2 + 1e-24)));
            localQ = (GramMatrix + RegParam^2*TVterm);
            localDistribution = 0*localDistribution;
            for i=1:nSignals
                localDistribution = localDistribution + weights(i)*localQ\Kernel{i}.'*Signal{i};
            end
            change = norm(localDistribution - prev);
            %Stop if result is stable
            if round(change,5) ==0
                break;
            end
        end
        %Compute the partial pseudinverse from the optimized result
        TVterm = RegMatrix.'*((RegMatrix./sqrt((RegMatrix*localDistribution).^2 + 1e-24)));
        Q = (GramMatrix + RegParam^2*TVterm);
        
    case 'huber'
        localDistribution = zeros(distDim,1);
        for j=1:500
            prev = localDistribution;
            %Compute pseudoinverse and unconst. distribution recursively
            HuberTerm = 1/(HuberParam^2)*(RegMatrix.'*(RegMatrix./sqrt((RegMatrix*localDistribution/HuberParam).^2 + 1)));
            localQ = (GramMatrix + RegParam^2*HuberTerm);
            localDistribution = 0*localDistribution;
            for i=1:nSignals
                localDistribution = localDistribution + weights(i)*localQ\Kernel{i}.'*Signal{i};
            end
            change = norm(localDistribution - prev);
            %Stop if result is stable
            if round(change,5) ==0
                break;
            end
        end
        %Compute the partial pseudinverse from the optimized result
        HuberTerm = 1/(HuberParam^2)*(RegMatrix.'*(RegMatrix./sqrt((RegMatrix*localDistribution/HuberParam).^2 + 1)));
        Q = (GramMatrix + RegParam^2*HuberTerm);
end

end