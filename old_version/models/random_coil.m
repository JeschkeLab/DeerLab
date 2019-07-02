function distr=random_coil(r0,par)
%
% Model library of DeerAnalysis2008: random_coil
%
% Random-coil model for an unfolded peptide/protein, end-to-end distance
% distribution is approximated by a Gaussian coil with proper mean
% distance, which is good for sufficiently large N
% N. C. Fitzkee, G. D. Rose, PNAS 2004, 101(34), 12497-12502
% equation (1) (has a typo, R has to be replaced by N) and value for R0 and
% default value for nu from figure caption Figure 4
%
% (c) G. Jeschke, 2008
%
% PARAMETERS
% name    symbol default lower bound upper bound
% par(1)  N      50      2              1000    number of residues between labels, including labeled residues
% par(2)  nu     0.602   0.33           1       scaling exponent

R0=0.198; % 1.98 Å per residue
N=par(1);
nu=par(2);

Rg=R0*N^nu;
R2=6*Rg^2; % mean square end-to-end distance from radius of gyration
c0=3/(2*pi*R2)^(3/2); % normalization prefactor
shell=4*pi*r0.^2; % spherical shell surface
garg=3*r0.^2/(2*R2); % argument for Gaussian distribution
gauss=exp(-garg);
distr=c0*shell.*gauss;
distr = distr/um(distr);
