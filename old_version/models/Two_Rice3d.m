function distr=Two_Rice3d(r0,par)
%
% Model for DEERAnalysis 2010: Two Rice3d
%
% Two Rice distributions (Rician distributions) with mean distances <nu1>,
% <nu2> and widths (standard deviations) sigma1, sigma2.
% (c) M. Spitzbarth, 2010
%
% See doi:10.1016/j.jmr.2010.10.005 for the formula.
%
% -----
% The mean distance ratio <nu2>/<nu1> is constrained by the lower 
% and upper ratio limits kmin, kmax. To disable use of distance 
% ratio constraint(s), you can 1) uncheck the desired kmin and/or
% kmax parameter checkbox or 2) delete the header's kmin and kmax
% parameter definition rows below.  
%
% H.C. Hyde, 2011
% -----
%
% PARAMETERS
% name    symbol default lower bound upper bound
% par(1)  <nu1>    2.5     1.0        10         mean distance
% par(2)  sigma1   0.4     0.1        5          standard deviation
% par(3)  <nu2>    4.0     1.0        10         mean distance
% par(4)  sigma2   0.4     0.1        5          standard deviation
% par(5)  p1       0.5     0          1          fraction of pairs at 1st distance
% par(6)  kmin     1.00    0          Inf        <nu2>/<nu1> minimum ratio
% par(7)  kmax     2.00    0          Inf        <nu2>/<nu1> maximum ratio


nu1=par(1);
sigma1=par(2);
nu2=par(3);
sigma2=par(4);
p1=par(5);

distr1=r0/nu1/(sigma1*sqrt(2*pi)).*( ...
    exp(-.5*((r0-nu1)/(sqrt(2)*sigma1)).^2) ...
    - exp(-.5*((r0+nu1)/(sqrt(2)*sigma1)).^2));
% The next line effectively multiplies the above with the heaviside
% step function. The Rice distribution is zero for negative values.
distr1(distr1<0)=0;

distr2=r0/nu2/(sigma2*sqrt(2*pi)).*( ...
    exp(-.5*((r0-nu2)/(sqrt(2)*sigma2)).^2) ...
    - exp(-.5*((r0+nu2)/(sqrt(2)*sigma2)).^2));
% The next line effectively multiplies the above with the heaviside
% step function. The Rice distribution is zero for negative values.
distr2(distr2<0)=0;

distr=p1*distr1+(1-p1)*distr2;
distr = distr/sum(distr);
