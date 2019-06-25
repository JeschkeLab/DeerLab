function distr=Gaussian_hom(r0,par)
%
% Model library of DeerAnalysis2006: Gaussian_hom
%
% single Gaussian peak with mean distance <r> and width (standard
% deviation) s(r) and homogeneous background (three-dimensional) with
% concentration c in arbitrary units
%
% (c) G. Jeschke, 2006
%
% PARAMETERS
% name    symbol default lower bound upper bound
% par(1)  <r>    3.5     1.5         10
% par(2)  s(r)   0.5     0.05        5
% par(3)  c      0.1     0           10

gauss0=(r0-par(1)*ones(size(r0)))/(sqrt(2)*par(2));
distr=exp(-gauss0.^2);
r2=r0.^2;
fac=sum(distr)/1e3; % normalization factor for hom. background rel to Gaussian
hom=par(3)*fac*r2;
distr=distr+hom;
distr = distr/sum(distr);
