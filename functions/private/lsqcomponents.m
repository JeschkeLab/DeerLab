%
% LSQCOMPONENTS Construct input arguments for NNLS solvers
%
%   [Q,KtV] = LSQCOMPONENTS(V,r,K,L,alpha,'type')
%   Computes the components required by non-negative least-squares (NNLS)
%   solvers (Q) and (KtV) for fitting a regularization model. The function
%   requires the signal (V), dipolar kernel (K), regularization matrix (L)
%   and regularization parameter (alpha). The type of regularization
%   functional is determined by the 'type' string argument.
%
%   [Q,KtV] = LSQCOMPONENTS(V,r,K,L,alpha,'type',eta)
%   The Huber parameter can be specified by passing as the (eta) argument.
%
%   [Q,KtV,w] = LSQCOMPONENTS({V1,V2,...},r,{K1,K2,...},L,alpha,'type',eta)
%   Passing multiple signals/kernels constructs the LSQ components (Q) and (KtS)
%   as required for global fit of the regularization functionals. The global fit
%   weights are automatically computed according to their contribution
%   to ill-posedness. The calculated weights (w) can be requested as an
%   additional output argument.
%
%   [Q,KtV] = LSQCOMPONENTS(V,r,K,L,alpha,'type',eta,w)
%   The global fit weights (w) can be manually passed to avoid computing them
%   automatically.
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function [Q,KtV,weights] = lsqcomponents(V,r,K,L,alpha,RegType,HuberParam,weights)

% Ensure that signals and kernel are in a cell array
if ~iscell(V)
    V = {V};
end
if ~iscell(K)
    K = {K};
end
% Provide defaults for Huber parameter and weights
if nargin<7 || isempty(HuberParam)
   HuberParam = 1.35; 
end
if nargin<8 
    weights = [];
end

% Prepare
nSignals = length(V);
distDim = length(L);

% Get weights of different signals for global fitting
if isempty(weights)
    weights = globalweights(V);
end

% Compute the terms depending on different signals
KtV = zeros(distDim,1);
KtK = zeros(distDim,distDim);
for i = 1:nSignals
    KtV = KtV + weights(i)*K{i}.'*V{i};
    KtK = KtK + weights(i)*K{i}.'*K{i};
end

% Compute then the LSQ components needed by NNLS optimizers
switch lower(RegType)
    
    case 'tikhonov'
        Q = KtK + alpha^2*(L.'*L);
        
    case 'tv'
        P = zeros(distDim,1);
        maxIter = 500;
        TVterm = @(p)L.'*((L./sqrt((L*p).^2 + eps)));
        for j = 1:maxIter
            Pprev = P;
            % Compute pseudoinverse and unconst. distribution recursively 
            Q_ = KtK + alpha^2*TVterm(P);
            P = zeros(distDim,1);
            for i = 1:nSignals
                P = P + weights(i)*Q_\K{i}.'*V{i};
            end
            % Normalize distribution by its integral to stabilize convergence
            P = P/sum(abs(P))/mean(diff(r));
            % Monitor largest changes in the distribution
            change = max(abs(P - Pprev));    
                        
            % Stop if result is stable
            if change < 1e-1
                break;
            end
        end
        % Compute the partial pseudoinverse from the optimized result
        Q = KtK + alpha^2*TVterm(P);
        
    case 'huber'
        P = zeros(distDim,1);
        maxIter = 500;
        HuberTerm = @(p) 1/(HuberParam^2)*(L.'*(L./sqrt((L*p/HuberParam).^2 + 1)));
        for j = 1:maxIter
            Pprev = P;
            % Compute pseudoinverse and unconst. distribution recursively
            Q_ = KtK + alpha^2*HuberTerm(P);
            P = zeros(distDim,1);
            for i = 1:nSignals
                P = P + weights(i)*Q_\K{i}.'*V{i};
            end
            % Normalize distribution by its integral to stabilize convergence
            P = P/sum(abs(P))/mean(diff(r));
            % Monitor largest changes in the distribution
            change = max(abs(P - Pprev));
            % Stop if result is stable
            if change < 1e-2
                break;
            end
        end
        % Compute the partial pseudoinverse from the optimized result
        Q = KtK + alpha^2*HuberTerm(P);
end

end
