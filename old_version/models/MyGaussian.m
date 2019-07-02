function distr = MyGaussian(r0,par)
%
% Model library of DeerAnalysis2006: Gaussian
%
% single Gaussian peak with mean distance <r> and width (standard
% deviation) s(r)
% (c) G. Jeschke, 2006,2019
%
% PARAMETERS
% name    symbol default lower bound upper bound
% par(1)  <r>    3.5     1.0         10         mean distance
% par(2)  s(r)   0.5     0.2         0.5        standard deviation

gauss0=(r0-par(1)*ones(size(r0)))/(sqrt(2)*par(2));
distr=exp(-gauss0.^2);
