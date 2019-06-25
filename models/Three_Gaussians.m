function distr = Three_Gaussians(r0,par)
%
% Model library of DeerAnalysis2018: Three_Gaussians
%
% Two Gaussian peaks with mean distances <r1>, <r2>, <r3> 
% and widths (standard deviation) s(r1), s(r2), s(r3)
% the peaks have relative integral intensities
% p1, p2, and p3 = 1 - p1 - p2.
%
% G. Jeschke, 2018
%
% -----
%
% PARAMETERS
% name    symbol default lower bound upper bound
% par(1)  <r1>   2.5     1.5         20         1st mean distance
% par(2)  s(r1)  0.5     0.05        5          std. dev. of 1st distance
% par(3)  <r2>   3.5     1.5         20         2nd mean distance
% par(4)  s(r2)  0.5     0.05        5          std. dev. of 2nd distance
% par(5)  <r3>   5.0     1.5         20         3rd mean distance
% par(6)  s(r3)  0.5     0.05        5          std. dev. of 3rd distance
% par(7)  p1     0.5     0           1          fraction of pairs at 1st
%                                               distance
% par(8)  p2     0.3     0           1          fraction of pairs at 2nd
%                                               distance

gauss1=(r0-par(1)*ones(size(r0)))/(sqrt(2)*par(2));
gauss1=exp(-gauss1.^2);
intg1=sum(gauss1);
gauss2=(r0-par(3)*ones(size(r0)))/(sqrt(2)*par(4));
gauss2=exp(-gauss2.^2);
gauss2=gauss2*intg1/sum(gauss2);
gauss3=(r0-par(5)*ones(size(r0)))/(sqrt(2)*par(6));
gauss3=exp(-gauss3.^2);
gauss3=gauss3*intg1/sum(gauss3);
distr=par(7)*gauss1+par(8)*gauss2+(1-par(7)-par(8))*gauss3;
distr = distr/sum(distr);

