function [idx_Lcorner,idx_AIC,idx_GCV,rho,eta,alpha] = l_curve_mod(K,L,S,noise)
%L_CURVE_MOD Solve Tikhonov for a range of regularization parameters and
% determine the optimal ones based on several criteria.

% Set defaults
if ~exist('noise','var'), noise = 0; end
unconstrainedTikhonov = true;

% Preparations
L = full(L);
nt = numel(S);
KtK = K.'*K;
LtL = L.'*L;
KtS = K.'*S;

% Get alpha range
alpha = get_regparamrange(K,L,noise);

% Calculate all metrics over alpha range
%-------------------------------------------------------------
for a = numel(alpha):-1:1
  Q = KtK + alpha(a)^2*LtL;
  if unconstrainedTikhonov
    P = Q\KtS; % unconstrained Tikhonov solution
  else
    P = fnnls(Q,KtS); % non-negative Tikhonov solution
  end
  Serr = K*P - S; % time-domain fit residuals
  H = K*(Q\K.'); % influence/projection/hat matrix
  rho(a) = log(norm(Serr)); % time-domain fit error
  eta(a) = log(norm(L*P)); % distance-domain roughness
  AIC_metric(a) = nt*log(norm(Serr)^2/nt) + 2*trace(H);
  GCV_metric(a) = norm(Serr)^2/(1-trace(H)/nt)^2;
end
rescaleL = @(x) (x-min(x))/(max(x)-min(x));
Lcorner_metric = rescaleL(eta).^2 + rescaleL(rho).^2;

% Find alpha values at AIC, GCV, and Lcorner metric minima
%-------------------------------------------------------------
[~,idx_AIC] = min(AIC_metric);
[~,idx_GCV] = min(GCV_metric);
[~,idx_Lcorner] = min(Lcorner_metric);

% Print results
%-------------------------------------------------------------
printResults = false;
if printResults
  fprintf('Regularization parameter search:\n');
  fprintf('  range:    alpha = %g to %g  log10(alpha) = %g to %g  (%d values)\n',...
    alpha(1),alpha(end),log10(alpha(1)),log10(alpha(end)),numel(alpha));
  fprintf('  AIC:      alpha = %g  log10(alpha) = %g  (idx %d)\n',...
    alpha(idx_AIC),log10(alpha(idx_AIC)),idx_AIC);
  fprintf('  GCV:      alpha = %g  log10(alpha) = %g  (idx %d)\n',...
    alpha(idx_GCV),log10(alpha(idx_GCV)),idx_GCV);
  fprintf('  L-corner: alpha = %g  log10(alpha) = %g  (idx %d)\n',...
    alpha(idx_Lcorner),log10(alpha(idx_Lcorner)),idx_Lcorner);
end
