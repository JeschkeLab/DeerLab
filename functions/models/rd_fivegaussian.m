%
% RD_FIVEGAUSSIAN Sum of four Gaussian distributions parametric model
%
%   info = RD_FIVEGAUSSIAN
%   Returns an (info) structure containing the specifics of the model.
%
%   P = RD_FIVEGAUSSIAN(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the paramteres array (param). The required parameters can also be found 
%   in the (info) structure.
%
% PARAMETERS
% name      symbol default lower bound upper bound
% --------------------------------------------------------------------------
% param(1)  <r1>   2.5     1.5         20         1st mean distance
% param(2)  w1     0.5     0.2         5          FWHM of 1st distance
% param(3)  <r2>   3.0     1.5         20         2nd mean distance
% param(4)  w2     0.5     0.2         5          FWHM of 2nd distance
% param(5)  <r3>   3.5     1.5         20         3rd mean distance
% param(6)  w3     0.5     0.2         5          FWHM of 3rd distance
% param(7)  <r4>   4.5     1.5         20         4th mean distance
% param(8)  w4     0.5     0.2         5          FWHM of 4th distance
% param(9)  <r5>   5.0     1.5         20         5th mean distance
% param(10) w5     0.5     0.2         5          FWHM of 5th distance
% param(11) p1     0.2     0           1          Fraction of pairs at 1st distance
% param(12) p2     0.2     0           1          Fraction of pairs at 2nd distance
% param(13) p3     0.2     0           1          Fraction of pairs at 3rd distance
% param(14) p4     0.2     0           1          Fraction of pairs at 4th distance
%--------------------------------------------------------------------------
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function output = rd_fivegaussian(r,param)

nParam = 14;

if nargin~=0 && nargin~=2
    error('Model requires two input arguments.')
end

if nargin==0
    %If no inputs given, return info about the parametric model
    info.model  = 'Five-Gaussian distribution';
    info.nparam  = nParam;
    info.parameters(1).name = 'Mean distance <r1> of 1st Gaussian';
    info.parameters(1).range = [1 20];
    info.parameters(1).default = 2.5;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = 'FWHM w1 of 1st Gaussian';
    info.parameters(2).range = [0.2 5];
    info.parameters(2).default = 0.5;
    info.parameters(2).units = 'nm';
    
    info.parameters(3).name = 'Mean distance <r2>  of 2nd Gaussian';
    info.parameters(3).range = [1 20];
    info.parameters(3).default = 3;
    info.parameters(3).units = 'nm';
    
    info.parameters(4).name = 'FWHM w2 of 2nd Gaussian';
    info.parameters(4).range = [0.2 5];
    info.parameters(4).default = 0.5;
    info.parameters(4).units = 'nm';
    
    info.parameters(5).name = 'Mean distance <r3> of 3rd Gaussian';
    info.parameters(5).range = [1 20];
    info.parameters(5).default = 3.5;
    info.parameters(5).units = 'nm';
    
    info.parameters(6).name = 'FWHM w3 of 3rd Gaussian';
    info.parameters(6).range = [0.2 5];
    info.parameters(6).default = 0.5;
    info.parameters(6).units = 'nm';
    
    info.parameters(7).name = 'Mean distance <r4> of 4th Gaussian';
    info.parameters(7).range = [1 20];
    info.parameters(7).default = 4;
    info.parameters(7).units = 'nm';
    
    info.parameters(8).name = 'FWHM w4 of 4th Gaussian';
    info.parameters(8).range = [0.2 5];
    info.parameters(8).default = 0.5;
    info.parameters(8).units = 'nm';
    
    info.parameters(9).name = 'Mean distance <r4> of 4th Gaussian';
    info.parameters(9).range = [1 20];
    info.parameters(9).default = 5;
    info.parameters(9).units = 'nm';
    
    info.parameters(10).name = 'FWHM w4 of 4th Gaussian';
    info.parameters(10).range = [0.2 5];
    info.parameters(10).default = 0.5;
    info.parameters(10).units = 'nm';
    
    info.parameters(11).name = 'Relative amplitude A1 of 1st Gaussian';
    info.parameters(11).range = [0 1];
    info.parameters(11).default = 0.2;
    
    info.parameters(12).name = 'Relative amplitude A2 of 2nd Gaussian';
    info.parameters(12).range = [0 1];
    info.parameters(12).default = 0.2;
    
    info.parameters(13).name = 'Relative amplitude A3 of 3rd Gaussian';
    info.parameters(13).range = [0 1];
    info.parameters(13).default = 0.2;
       
    info.parameters(14).name = 'Relative amplitude A4 of 4th Gaussian';
    info.parameters(14).range = [0 1];
    info.parameters(14).default = 0.2;
    
    output = info;
    return
end

if nargin~=2
    error('Model requires two input arguments.')
end

% Compute the model distance distribution
Lam1 = (param(2)/sqrt(2*log(2)));
Lam2 = (param(4)/sqrt(2*log(2)));
Lam3 = (param(6)/sqrt(2*log(2)));
Lam4 = (param(8)/sqrt(2*log(2)));
Lam5 = (param(10)/sqrt(2*log(2)));
Gaussian1 = sqrt(2/pi)*1/Lam1*exp(-2*((r(:) - param(1))/Lam1).^2);
Gaussian2 = sqrt(2/pi)*1/Lam2*exp(-2*((r(:) - param(3))/Lam2).^2);
Gaussian3 = sqrt(2/pi)*1/Lam3*exp(-2*((r(:) - param(5))/Lam3).^2);
Gaussian4 = sqrt(2/pi)*1/Lam4*exp(-2*((r(:) - param(7))/Lam4).^2);
Gaussian5 = sqrt(2/pi)*1/Lam5*exp(-2*((r(:) - param(9))/Lam5).^2);
P = param(11)*Gaussian1 + param(12)*Gaussian2 + param(13)*Gaussian3 ...
     + param(14)*Gaussian4 +max(1 - param(10) - param(12) - param(13) - param(14),0)*Gaussian5;

% Normalize
dr = r(2)-r(1);
P = P/sum(P)/dr;

output = P;

return