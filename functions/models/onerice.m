function output = onerice(r,param)

% PARAMETERS
% name    symbol default lower bound upper bound
% par(1)  <nu>    3.5     1.0         10          mean distance
% par(2)  sigma   0.7     0.1          5          standard deviation


nParam = 2;

if nargin==0
    %If no inputs given, return info about the parametric model
    info.Model  = 'Single Rice/Rician distribution';
    info.Equation  = ['x/',char(963),'²*exp((r-',char(957),')²/(2',char(963),'²))*Bessel(r*',char(957),'/',char(963),'²)'];
    info.nParam  = nParam;
    info.parameters(1).name = ['Mean distance ',char(957)];
    info.parameters(1).range = [1 10];
    info.parameters(1).default = 3.5;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = ['Standard deviation ',char(963)];
    info.parameters(2).range = [0.1 5];
    info.parameters(2).default = 0.7;
    info.parameters(2).units = 'nm';
    
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
    Distribution = (r./sqscale).*exp(-1/2*(r.^2 + nu.^2)./sqscale).*besseli_(0,r.*nu./sqscale);
    %The Rice distribution is zero for negative values.
    Distribution(Distribution<0)=0;
    
    if ~iscolumn(Distribution)
        Distribution = Distribution';
    end
    if ~all(Distribution==0)
        Distribution = Distribution/sum(Distribution);
    end
    output = Distribution;
else
    
    %Else, the user has given wrong number of inputs
    error('Model requires two input arguments.')
end


return