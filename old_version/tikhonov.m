function P_alpha = tikhonov(K,L,S,alpha)
% TIKHONOV   Tikhonov regularization
%
%  P_alpha = tikhonov(K,L,S,alpha)
%
% Computes the unconstrained solution
%
%   P_alpha = argmin(||K*P-S||^2 + alpha^2*||L*P||^2)
%
% Input:
%    K      kernel matrix, N x N
%    L      regularization matrix, M x N
%    S      time-domain data vector (form factor), N x 1
%    alpha  regularization parameter (scalar or vector)
%
% Output:
%    P_alpha  distance-domain solution vector, N x 1 (if a one alpha is given)
%             matrix of solution vectors, N x Q (if Q alpha values are given)

KtK = K.'*K;
LtL = L.'*L;
KtS = K.'*S;
for a = numel(alpha):-1:1
  Q = KtK + alpha(a)^2*LtL;
  P_alpha(:,a) = Q\KtS;
end
