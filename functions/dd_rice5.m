%
% DD_RICE5 Sum of five 3D-Rice distributions parametric model
%
%   info = DD_RICE5
%   Returns an (info) structure containing the specifics of the model.
%
%   P = DD_RICE5(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to
%   the paramteres array (param). The required parameters can also be found
%   in the (info) structure.
%
% PARAMETERS
% name      symbol default lower bound upper bound
% --------------------------------------------------------------------------
% param(1)   nu1      2.5     1.0        10         center of 1st component
% param(2)   sigma1   0.7     0.1        5          spread of 1st component
% param(3)   nu2      3.5     1.0        10         center of 2nd component
% param(4)   sigma2   0.7     0.1        5          spread of 2nd component
% param(5)   nu3      4.0     1.0        10         center of 3rd component
% param(6)   sigma3   0.7     0.1        5          spread of 3rd component
% param(7)   nu4      5.0     1.0        10         center of 4th component
% param(8)   sigma4   0.7     0.1        5          spread of 4th component
% param(9)   nu5      5.5     1.0        10         center of 5th component
% param(10)  sigma5   0.7     0.1        5          spread of 5th component
% param(11)  p1       0.2     0          1          amplitude of 1st component
% param(12)  p2       0.2     0          1          amplitude of 2nd component
% param(13)  p3       0.2     0          1          amplitude of 3rd component
% param(14)  p3       0.2     0          1          amplitude of 4th component
% --------------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function output = dd_rice5(r,param)


nParam = 14;

if nargin~=0 && nargin~=2
    error('Model requires two input arguments.')
end

if nargin==0
    % If no inputs given, return info about the parametric model
    info.model  = 'Five 3D-Rice distributions';
    info.nparam  = nParam;
    info.parameters(1).name = ['Center ',char(957),'1 1st component'];
    info.parameters(1).range = [1 10];
    info.parameters(1).default = 2.5;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = ['Spread ',char(963),'1 1st component'];
    info.parameters(2).range = [0.1 5];
    info.parameters(2).default = 0.7;
    info.parameters(2).units = 'nm';
    
    info.parameters(3).name = ['Center ',char(957),'2 2nd component'];
    info.parameters(3).range = [1 10];
    info.parameters(3).default = 3.5;
    info.parameters(3).units = 'nm';
    
    info.parameters(4).name = ['Spread ',char(963),'2 2nd component'];
    info.parameters(4).range = [0.1 5];
    info.parameters(4).default = 0.7;
    info.parameters(4).units = 'nm';
    
    info.parameters(5).name = ['Center ',char(957),'3 3rd component'];
    info.parameters(5).range = [1 10];
    info.parameters(5).default = 4.0;
    info.parameters(5).units = 'nm';
    
    info.parameters(6).name = ['Spread ',char(963),'3 3rd component'];
    info.parameters(6).range = [0.1 5];
    info.parameters(6).default = 0.7;
    info.parameters(6).units = 'nm';
    
    info.parameters(7).name = ['Center ',char(957),'4 4th component'];
    info.parameters(7).range = [1 10];
    info.parameters(7).default = 5;
    info.parameters(7).units = 'nm';
    
    info.parameters(8).name = ['Spread ',char(963),'4 4th component'];
    info.parameters(8).range = [0.1 5];
    info.parameters(8).default = 0.7;
    info.parameters(8).units = 'nm';
    
    info.parameters(9).name = ['Center ',char(957),'5 5th component'];
    info.parameters(9).range = [1 10];
    info.parameters(9).default = 5.5;
    info.parameters(9).units = 'nm';
    
    info.parameters(10).name = ['Spread ',char(963),'5 5th component'];
    info.parameters(10).range = [0.1 5];
    info.parameters(10).default = 0.7;
    info.parameters(10).units = 'nm';
    
    info.parameters(11).name = 'Relative amplitude A1 1st component';
    info.parameters(11).range = [0 1];
    info.parameters(11).default = 0.2;
    info.parameters(11).units = '';
    
    info.parameters(12).name = 'Relative amplitude A2 2nd component';
    info.parameters(12).range = [0 1];
    info.parameters(12).default = 0.2;
    info.parameters(12).units = '';
     
    info.parameters(13).name = 'Relative amplitude A3 3rd component';
    info.parameters(13).range = [0 1];
    info.parameters(13).default = 0.2;
    info.parameters(13).units = '';
     
    info.parameters(14).name = 'Relative amplitude A4 4th component';
    info.parameters(14).range = [0 1];
    info.parameters(14).default = 0.2;
    info.parameters(14).units = '';
        
    output = info;
    return;
end

% If user passes them, check that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

% Parse input
validateattributes(r,{'numeric'},{'nonnegative','increasing','nonempty'},mfilename,'r')

% Compute non-central chi distribution with 3 degrees of freedom (a 3D Rician)
nu = param([1 3 5 7 9]);
sig = param([2 4 6 8 10]);
a = param([10 11 12 13]);
a(5) = max(1-sum(a),0);

P = multirice3d(r,nu,sig,a);

output = P;

return
