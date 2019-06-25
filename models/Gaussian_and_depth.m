function [sim,distr]=Gaussian_and_depth(r0,t,par)
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
% par(2)  s(r)   0.5     0.02        5          standard deviation
% par(3)  Delta  0.5     0.001       1.0        modulation depth

gauss0=(r0-par(1)*ones(size(r0)))/(sqrt(2)*par(2));
distr=exp(-gauss0.^2);

[sim,sc]=deer_sim(r0,distr,t,20);
sim=sim-0.99*ones(size(sim));
sim=sim*par(3)/0.01;
sim=sim+(1-par(3))*ones(size(sim));