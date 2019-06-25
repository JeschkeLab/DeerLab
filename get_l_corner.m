function poi=get_l_corner(rho,eta)
%
% Determines optimum regularization parameter from discrete L curve
% (maximum curvature)
%
% Input:
% rho       log(square of residual norm)
% eta       log(square of (semi)norm of 2nd derivative of solution)



etamin=min(eta);
rhomin=min(rho);
etamax=max(eta);
rhomax=max(rho);
eta1=(eta-etamin)/(etamax-etamin);
rho1=(rho-rhomin)/(rhomax-rhomin);

[mi,poi]=min(eta1.^2+rho1.^2);

