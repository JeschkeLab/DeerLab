function output = tworice(r,param)

% PARAMETERS
% name    symbol default lower bound upper bound
% par(1)  <nu1>    2.5     1.0        10         mean distance
% par(2)  sigma1   0.4     0.1        5          standard deviation
% par(3)  <nu2>    4.0     1.0        10         mean distance
% par(4)  sigma2   0.4     0.1        5          standard deviation
% par(5)  p1       0.5     0          1          fraction of pairs at 1st distance
% par(6)  kmin     1.00    0          Inf        <nu2>/<nu1> minimum ratio
% par(7)  kmax     2.00    0          Inf        <nu2>/<nu1> maximum ratio


nParam = 5;

if nargin==0
    %If no inputs given, return info about the parametric model
    info.Model  = 'Two Rice/Rician distributions';
    info.Equation  = ['A*r/',char(963),'1²*exp((r-',char(957),'1)²/(2',char(963),'1²))*Bessel(r*',char(957),'1/',char(963),'1²)'...
        '+ (1-A)*r/',char(963),'2²*exp((r-',char(957),'2)²/(2',char(963),'2²))*Bessel(r*',char(957),'2/',char(963),'2²)' ];
    info.nParam  = nParam;
    info.parameters(1).name = ['Mean distance ',char(957),'1 1st Rician'];
    info.parameters(1).range = [1 10];
    info.parameters(1).default = 2.5;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = ['Standard deviation ',char(963),'1 1st Rician'];
    info.parameters(2).range = [0.1 5];
    info.parameters(2).default = 0.7;
    info.parameters(2).units = 'nm';
    
    info.parameters(3).name = ['Mean distance ',char(957),'2 2nd Rician'];
    info.parameters(3).range = [1 10];
    info.parameters(3).default = 4.5;
    info.parameters(3).units = 'nm';
    
    info.parameters(4).name = ['Standard deviation ',char(963),'2 2nd Rician'];
    info.parameters(4).range = [0.1 5];
    info.parameters(4).default = 0.7;
    info.parameters(4).units = 'nm';
    
    info.parameters(5).name = 'Relative amplitude A 1st Rician';
    info.parameters(5).range = [0 1];
    info.parameters(5).default = 0.5;
    
    output = info;
    
elseif nargin == 2
    
    %If user passes them, check that the number of parameters matches the model
    if length(param)~=nParam
        error('The number of input parameters does not match the number of model parameters.')
    end
    
    nu=param(1);
    sqscale=param(2).^2;
    %Compute rician/rice distribution using the zeroth order modified Bessel function of
    %the first kind
    Rician1 = (r./sqscale).*exp(-1/2*(r.^2 + nu.^2)./sqscale).*besseli_(0,r.*nu./sqscale);
    %The Rice distribution is zero for negative values.
    Rician1(Rician1<0)=0;
    
    nu=param(3);
    sqscale=param(4).^2;
    %Compute rician/rice distribution using the zeroth order modified Bessel function of
    %the first kind
    Rician2 = (r./sqscale).*exp(-1/2*(r.^2 + nu.^2)./sqscale).*besseli_(0,r.*nu./sqscale);
    %The Rice distribution is zero for negative values.
    Rician2(Rician2<0)=0;
    
    %Construct distance distribution
    Distribution = param(5)*Rician1 + max(1-param(5),0)*Rician2;
    
    if ~iscolumn(Distribution)
        Distribution = Distribution';
    end
    if ~all(Distribution==0)
        Distribution = Distribution/sum(Distribution)/mean(diff(r));
    end
    output = Distribution;
else
    
    %Else, the user has given wrong number of inputs
    error('Model requires two input arguments.')
end


return