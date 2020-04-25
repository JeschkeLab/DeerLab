%
% LSQCOMPONENTS Construct input arguments for NNLS solvers
%
%   [KtKreg,KtV] = LSQCOMPONENTS(V,K,r,L,alpha,regtype)
%   [KtKreg,KtV] = LSQCOMPONENTS(V,K,r,L,alpha,regtype,eta)
%   [KtKreg,KtV] = LSQCOMPONENTS({V1,V2,...},{K1,K2,...},r,L,alpha,regtype,eta,w)
%
%   Computes the components required by non-negative least-squares (NNLS)
%   solvers (KtKreg) and (KtV) for fitting a regularization model. The function
%   requires the signal (V), the dipolar kernel (K), the distance axis (r),
%   the regularization matrix (L), and the regularization parameter (alpha).
%
%   The type of regularization functional is determined by (regtype).
%   The Huber parameter can be specified by passing as the (eta) argument.
%
%   Passing multiple signals/kernels constructs the LSQ components (KtKreg) and
%   (KtV) as required for global fit of the regularization functionals. In this
%   case, the weights need to passed in (w).
%
% Inputs:
%   V       signal
%   K       kernel matrix (us and nm)
%   r       distance axis (nm)
%   L       regularization operator matrix (see regoperator)
%   alpha   regularization parameter
%   regtype regularzion type: 'Tikhonovov', 'Huber', 'TV'
%   eta     parameter for Huber regularization (optional, default 1.35)
%   w       global weights (required only if V has multiple signals)
%
% Output:
%   KtKreg  K.'*K plus regularization term (weighted sum for multiple signals)
%   KtV     K.*V (weighted sum for multiple signals)
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function [KtKreg,KtV] = lsqcomponents(V,K,r,L,alpha,RegType,HuberParam,weights)

% Ensure that signals and kernel are in a cell array
if ~iscell(V)
    V = {V};
end
if ~iscell(K)
    K = {K};
end

% Prepare
nSignals = numel(V);
nr = size(L,2);

% Provide defaults for Huber parameter and weights
if nargin<7 || isempty(HuberParam)
   HuberParam = 1.35; 
end
if nargin<8 
    if nSignals==1
        weights = 1;
    else
        error('For multiple signals, weights need to be provided.');
    end
end

% Compute the terms as a weighted sum over components
KtV = zeros(nr,1);
KtK = zeros(nr,nr);
for i = 1:nSignals
    KtV = KtV + weights(i)*K{i}.'*V{i};
    KtK = KtK + weights(i)*K{i}.'*K{i};
end

% Compute then the LSQ components needed by NNLS optimizers
switch lower(RegType)
    
    case 'tikhonov'
        regterm = L.'*L;
        
    case 'tv'
        maxIter = 500;
        changeThreshold = 1e-1;
        TVFcn = @(p)L.'*((L./sqrt((L*p).^2 + eps)));
        regterm = optimizeregterm(TVFcn);
        
    case 'huber'
        maxIter = 500;
        changeThreshold = 1e-2;
        HuberFcn = @(p) 1/(HuberParam^2)*(L.'*(L./sqrt((L*p/HuberParam).^2 + 1)));
        regterm = optimizeregterm(HuberFcn);
end

KtKreg = KtK + alpha^2*regterm;

    function regterm = optimizeregterm(fun)
        P = zeros(nr,1);
        for j = 1:maxIter
            Pprev = P;
            % Compute pseudoinverse and unconst. distribution recursively
            KtKreg_ = KtK + alpha^2*fun(P);
            P = zeros(nr,1);
            for ii = 1:numel(V)
                P = P + weights(ii)*KtKreg_\K{ii}.'*V{ii};
            end
            % Normalize distribution by its integral to stabilize convergence
            P = P/sum(abs(P))/mean(diff(r));
            % Monitor largest changes in the distribution
            change = max(abs(P - Pprev));
            % Stop if result is stable
            if change < changeThreshold
                break;
            end
        end
        % Compute the regularization term from the optimized result
        regterm = fun(P);
    end


end

