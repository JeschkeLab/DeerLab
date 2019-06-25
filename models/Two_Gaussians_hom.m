function distr=Two_Gaussians_hom(r0,par)
%
% Model library of DeerAnalysis2006: Two_Gaussians_hom
%
% Two Gaussian peaks with mean distances <r1>, <r2> 
% and widths (standard deviation) s(r1), s(r2)
% the peaks have relative integral intensities
% p1 and p2=1-p1.
% A homogeneous background (three-dimensional) with
% concentration c in arbitrary units is added
%
% (c) G. Jeschke, 2006
%
% -----
% The mean distance ratio <r2>/<r1> is constrained by the lower 
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
% par(1)  <r1>   2.5     1.5         10
% par(2)  s(r1)  0.5     0.05        5
% par(3)  <r2>   3.5     1.5         10
% par(4)  s(r2)  0.5     0.05        5
% par(5)  p1     0.5     0           1
% par(6)  c      0.1     0           10
% par(7)  kmin   1.00    0           Inf        <r2>/<r1> minimum ratio   
% par(8)  kmax   2.00    0           Inf        <r2>/<r1> maximum ratio

gauss1=(r0-par(1)*ones(size(r0)))/(sqrt(2)*par(2));
gauss1=exp(-gauss1.^2);
intg1=sum(gauss1);
gauss2=(r0-par(3)*ones(size(r0)))/(sqrt(2)*par(4));
gauss2=exp(-gauss2.^2);
gauss2=gauss2*intg1/sum(gauss2);
distr=par(5)*gauss1+(1-par(5))*gauss2;
r2=r0.^2;
fac=sum(distr)/sum(r2); % normalization factor for hom. background rel to Gaussian
hom=par(6)*fac*r2;
distr=distr+hom;
distr = distr/sum(distr);
