function [rout,distr,rho,eta,reg_param,idx_corner,idx_AIC,idx_GCV] = get_Tikhonov_new(handles,reg_param)

persistent kdim K r t L

% Determine whether L curve should be calculated or not
calcLcurve = ~exist('reg_param','var') || length(reg_param)~=1;

tdip = handles.A_tdip;

if length(tdip) > 2048 % use standard kernel for very long data sets
    K = handles.Tikh_kernel;
    r = handles.Tikh_r;
    t = handles.Tikh_t;
    L = handles.Tikh_L;
    kdim = length(r);
    n1 = length(t);
    tf = linspace(0,max(tdip),n1);
    dt2 = tf(2)-tf(1);
    S = interp1(tdip,handles.A_dipevo,tf);
else
    if isempty(kdim) || kdim ~= length(tdip) % use previously stored kernel if dimension matches
        kdim = length(tdip);
        [K,t,r] = get_bas_Tikh(kdim);
        L = get_l(length(r),2); % differential regularization operator matrix
    end
    S = handles.A_dipevo;
    dt2 = tdip(2)-tdip(1);
end

dt = t(2)-t(1);
rout = r * (dt2/dt)^(1/3);

% % Code for testing
% figure(17); clf;
% plot(tdip,handles.A_dipevo,'k');
% hold on;
% plot(tf,ff,'r');
% fprintf(1,'New data size is (%i,%i).\n',size(ff));

S = S(:);
if calcLcurve
    % Compute L curve rho and eta (without nonnegativity constraint) and its corner
    [idx_corner,idx_AIC,idx_GCV,rho,eta,reg_param] = l_curve_mod(K,L,S,handles.fit_rms_value);
    alpha = reg_param(idx_corner);
else
    % Compute rho and eta for single regularization parameter
    Pfit = tikhonov(K,L,S,reg_param);
    rho = norm(K*Pfit-S);
    eta = norm(L*Pfit);
    alpha = reg_param(1);
    idx_corner = 1;
    idx_AIC = 1;
    idx_GCV = 1;
end

% Solve non-negativity-constrained Tikhonov problem
%---------------------------------------------------------------
% solver 1: active-set algorithm FNNLS (Matlab's lsqnonneg)
%    argmin(||C*x-d||^2) with C = [K;alpha*L] and d = [S;zeros]
% solver 2: active-set algorithm FNNLS (Matlab's lsqnonneg)
%    argmin(||S-K*distr||^2+alpha^2*||L*P||^2), with functional multiplied out
% solver 3: active-set algorithm FNNLS (own implementation)
%    argmin(||S-K*distr||^2+alpha^2*||L*P||^2), with functional multiplied out
% solver 4: block principal pivoting (own implementation)
%    argmin(||S-K*distr||^2+alpha^2*||L*P||^2), with functional multiplied out

nonnegSolver = 3;
switch nonnegSolver
  case 1
    C = [K;alpha*L];
    % options = optimset('Display','off','TolX',10*eps*max(size(C))*norm(kernel,1));
    options = optimset('Display','off','TolX',1000*eps);
    [m,~] = size(L);
    d = [S;zeros(m,1)];
    distr = lsqnonneg(C,d,options);
  case 2
    options = optimset('Display','off','TolX',1000*eps);
    Q = (K.'*K) + alpha^2*(L.'*L);
    distr = lsqnonneg(Q,K.'*S,options);
  case 3
    Q = (K.'*K) + alpha^2*(L.'*L);
    tol = 1e-8;
    distr = fnnls(Q,K.'*S,[],tol);
  case 4
    Q = (K.'*K) + alpha^2*(L.'*L);
    KtS = K.'*S;
    distr = nnls_bpp(Q,KtS,Q\KtS);
end
