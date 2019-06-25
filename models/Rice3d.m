function distr=Rice3d(r0,par)
%
% Model for DEERAnalysis 2010: Rice3d
%
% Rice distribution (Rician distribution) with mean distance <nu>
% and width (standard deviation) sigma.
% (c) M. Spitzbarth, 2010
%
% See doi:10.1016/j.jmr.2010.10.005 for the formula.
%
% PARAMETERS
% name    symbol default lower bound upper bound
% par(1)  <nu>    3.5     1.0         10          mean distance
% par(2)  sigma   0.7     0.1          5          standard deviation


nu=par(1);
sigma=par(2);

distr=r0/nu/(sigma*sqrt(2*pi)).*( ...
    exp(-.5*((r0-nu)/(sqrt(2)*sigma)).^2) ...
    - exp(-.5*((r0+nu)/(sqrt(2)*sigma)).^2));
% The next line effectively multiplies the above with the heaviside
% step function. The Rice distribution is zero for negative values.
distr(distr<0)=0;
distr = distr/sum(distr);