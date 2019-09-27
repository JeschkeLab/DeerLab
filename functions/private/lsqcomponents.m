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


function [Q,KtS,weights] = lsqcomponents(S,r,K,L,alpha,RegType,HuberParam,weights)

%Ensure that signals and kernel are in a cell array
if ~iscell(S)
    S = {S};
end
if ~iscell(K)
    K = {K};
end
%If Huber parameter not given, just use the default
if nargin<7 || isempty(HuberParam)
   HuberParam = 1.35; 
end
%If Huber parameter not given, just use the default
if nargin<8 
    weights = [];
end
%Prepare
nSignals = length(S);
distDim = length(L);
KtS = zeros(distDim,1);
GramMatrix = zeros(distDim,distDim);

%Get weights of different signals for global fitting
if isempty(weights)
weights = globalweights(S);
end
%Compute the terms depending on different signals
for i=1:nSignals
    KtS = KtS + weights(i)*K{i}.'*S{i};
    GramMatrix = GramMatrix + weights(i)*K{i}.'*K{i};
end

%Compute then the LSQ components needed by NNLS optimizers
switch lower(RegType)
    
    case 'tikhonov'
        Q = GramMatrix + alpha^2*(L.'*L);
        
    case 'tv'
        localP = zeros(distDim,1);
        for j=1:500
            prev = localP;
            %Compute pseudoinverse and unconst. distribution recursively
            TVterm = L.'*((L./sqrt((L*localP).^2 + eps)));
            localQ = (GramMatrix + alpha^2*TVterm);
            localP = 0*localP;
            for i=1:nSignals
                localP = localP + weights(i)*localQ\K{i}.'*S{i};
            end
            %Normalize distribution by its integral to stabilize convergence
            localP = localP/sum(abs(localP))/0.0099;
            %Monitor largest changes in the distribution
            change = max(localP - prev);    
                        
            %Stop if result is stable
            if change < 1e-1
                break;
            end
        end
        %Compute the partial pseudinverse from the optimized result
        TVterm = L.'*((L./sqrt((L*localP).^2 + 1e-24)));
        Q = (GramMatrix + alpha^2*TVterm);
        
    case 'huber'
        localP = zeros(distDim,1);
        for j=1:500
            prev = localP;
            %Compute pseudoinverse and unconst. distribution recursively
            HuberTerm = 1/(HuberParam^2)*(L.'*(L./sqrt((L*localP/HuberParam).^2 + 1)));
            localQ = (GramMatrix + alpha^2*HuberTerm);
            localP = 0*localP;
            for i=1:nSignals
                localP = localP + weights(i)*localQ\K{i}.'*S{i};
            end
            %Normalize distribution by its integral to stabilize convergence
            localP = localP/sum(abs(localP));
            %Monitor largest changes in the distribution
            change = max(localP - prev);
            %Stop if result is stable
            if change < 1e-2
                break;
            end
        end
        %Compute the partial pseudinverse from the optimized result
        HuberTerm = 1/(HuberParam^2)*(L.'*(L./sqrt((L*localP/HuberParam).^2 + 1)));
        Q = (GramMatrix + alpha^2*HuberTerm);
end

end