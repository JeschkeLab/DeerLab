function output = threegaussian(r,param)

% PARAMETERS
% name    symbol default lower bound upper bound
% par(1)  <r1>   2.5     1.5         20         1st mean distance
% par(2)  s(r1)  0.5     0.05        5          std. dev. of 1st distance
% par(3)  <r2>   3.5     1.5         20         2nd mean distance
% par(4)  s(r2)  0.5     0.05        5          std. dev. of 2nd distance
% par(5)  <r3>   5.0     1.5         20         3rd mean distance
% par(6)  s(r3)  0.5     0.05        5          std. dev. of 3rd distance
% par(7)  p1     0.5     0           1          fraction of pairs at 1st distance
% par(8)  p2     0.3     0           1          fraction of pairs at 2nd distance

nParam = 8;

if nargin==0
    %If no inputs given, return info about the parametric model
    info.Model  = 'Three-Gaussian distribution';
    info.Equation  = ['A1*exp(-(r-<r1>)²/(',char(963),'1*sqrt(2))²) + A2*exp(-(r-<r2>)²/(',char(963),'2*sqrt(2))²) + (1-A1-A2)*exp(-(r-<r3>)²/(',char(963),'3*sqrt(2))²)'];
    info.nParam  = nParam;
    info.parameters(1).name = 'Mean distance <r1> 1st Gaussian';
    info.parameters(1).range = [1 20];
    info.parameters(1).default = 2.5;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = ['Standard deviation ',char(963),'1 1st Gaussian'];
    info.parameters(2).range = [0.05 5];
    info.parameters(2).default = 0.5;
    info.parameters(2).units = 'nm';
    
    info.parameters(3).name = 'Mean distance <r2> 2nd Gaussian';
    info.parameters(3).range = [1 20];
    info.parameters(3).default = 3.5;
    info.parameters(3).units = 'nm';
    
    info.parameters(4).name = ['Standard deviation ',char(963),'2 2nd Gaussian'];
    info.parameters(4).range = [0.05 5];
    info.parameters(4).default = 0.5;
    info.parameters(4).units = 'nm';
    
    info.parameters(5).name = 'Mean distance <r3> 3rd Gaussian';
    info.parameters(5).range = [1 20];
    info.parameters(5).default = 3.5;
    info.parameters(5).units = 'nm';
    
    info.parameters(6).name = ['Standard deviation ',char(963),'3 3rd Gaussian'];
    info.parameters(6).range = [0.05 5];
    info.parameters(6).default = 0.5;
    info.parameters(6).units = 'nm';
    
    info.parameters(7).name = 'Relative amplitude A1 1st Gaussian';
    info.parameters(7).range = [0 1];
    info.parameters(7).default = 0.5;
    
    info.parameters(8).name = 'Relative amplitude A2 2nd Gaussian';
    info.parameters(8).range = [0 1];
    info.parameters(8).default = 0.5;
       
    output = info;
    
elseif nargin == 2
    
    %If user passes them, check that the number of parameters matches the model
    if length(param)~=nParam
        error('The number of input parameters does not match the number of model parameters.')
    end
    
    %If necessary inputs given, compute the model distance distribution
    Gaussian1 = exp(-((r-param(1))/(param(2))).^2);
    Gaussian2 = exp(-((r-param(3))/(param(4))).^2);
    Gaussian3 = exp(-((r-param(5))/(param(6))).^2);
    Distribution = param(7)*Gaussian1 + param(8)*Gaussian2 + max(1 - param(5) - param(8),0)*Gaussian3;
    if ~iscolumn(Distribution)
        Distribution = Distribution';
    end
    %Normalize
    Distribution = Distribution/sum(Distribution);
    output = Distribution;
else
    %Else, the user has given wrong number of inputs
    error('Model requires two input arguments.')
end

return