function functionHandle = regfunctional(Method,Signal,RegMatrix,Kernel,RegularizationParameter,HuberParameter)

if ~iscell(Signal)
    Signal = {Signal};
end
if ~iscell(Kernel)
    Kernel = {Kernel};
end
if nargin<6 || isempty(HuberParameter)
   HuberParameter = 1.35; 
end

%Get weights of different signals for global fitting
weights = globalweights(Signal);

switch Method
    case 'tv'
        functionHandle = @(Distribution)TVfunctional(Signal,Distribution,RegMatrix,Kernel,RegularizationParameter,weights);
    case 'tikhonov'
        functionHandle = @(Distribution)TikhonovFunctional(Signal,Distribution,RegMatrix,Kernel,RegularizationParameter,weights);
    case 'huber'
        functionHandle = @(Distribution)HuberFunctional(Signal,Distribution,RegMatrix,Kernel,RegularizationParameter,HuberParameter,weights);
end

end

%Tikhonov functional
%-----------------------------------------------------------
function [Functional,Gradient] = TikhonovFunctional(Signal,Distribution,L,Kernel,RegularizationParameter,weights)

Residual = 0;
ResidualGradient = 0;
for i=1:length(Signal)
    Residual = Residual + weights(i)*1/2*norm(Kernel{i}*Distribution - Signal{i})^2;
    ResidualGradient =  ResidualGradient + weights(i)*Kernel{i}.'*(Kernel{i}*Distribution - Signal{i});
end

Functional = Residual + RegularizationParameter^2*1/2*norm(L*Distribution)^2;
if nargout>1
    Gradient = ResidualGradient + RegularizationParameter^2*(L')*L*Distribution;
end

end

%Total variation functional
%-----------------------------------------------------------
function [Functional,Gradient] = TVfunctional(Signal,Distribution,L,Kernel,RegularizationParameter,weights)

Residual = 0;
ResidualGradient = 0;
for i=1:length(Signal)
    Residual = Residual + weights(i)*1/2*norm(Kernel{i}*Distribution - Signal{i})^2;
    ResidualGradient = ResidualGradient + weights(i)*Kernel{i}.'*(Kernel{i}*Distribution - Signal{i});
end
Functional = Residual + RegularizationParameter^2*sum(sqrt((L*Distribution).^2 + 1e-24));
if nargout>1
    TVGradient = (L)'*diag(1./(sqrt((L*Distribution).^2 + 1e-24)))*(L)*Distribution;
    Gradient = ResidualGradient + RegularizationParameter^2*TVGradient;
end

end

%Pseudo-Huber functional
%-----------------------------------------------------------
function [Functional,Gradient] = HuberFunctional(Signal,Distribution,L,Kernel,RegularizationParameter,HuberParameter,weights)

Residual = 0;
ResidualGradient = 0;
for i=1:length(Signal)
    Residual = Residual + weights(i)*1/2*norm(Kernel{i}*Distribution - Signal{i})^2;
    ResidualGradient = ResidualGradient + weights(i)*Kernel{i}.'*(Kernel{i}*Distribution - Signal{i});
end

HuberPenalty = (sqrt(1+(L*Distribution/HuberParameter).^2)-1);
HuberFunctional = sum(HuberPenalty);
%Get the Huber norm gradient
Functional = Residual + RegularizationParameter^2*HuberFunctional;
if nargout>1
    HuberGradient =  RegularizationParameter^2*(L/(HuberParameter))'*diag(1./sqrt((L*Distribution/HuberParameter).^2 + 1))*(L/(HuberParameter))*Distribution;
    Gradient = ResidualGradient + RegularizationParameter^2*HuberGradient;
end

end