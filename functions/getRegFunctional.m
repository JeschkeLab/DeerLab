function functionHandle = getRegFunctional(Method,Signal,RegMatrix,Kernel,RegularizationParameter,HuberParameter)

if nargin<6
    HuberParameter = 1.35;
end

switch Method
    case 'tv'
        functionHandle = @(Distribution)TVfunctional(Signal,Distribution,RegMatrix,Kernel,RegularizationParameter);
    case 'tikhonov'
        functionHandle = @(Distribution)TikhonovFunctional(Signal,Distribution,RegMatrix,Kernel,RegularizationParameter);
    case 'huber'
        functionHandle = @(Distribution)HuberFunctional(Signal,Distribution,RegMatrix,Kernel,RegularizationParameter,HuberParameter);
end

end

function [Functional,Gradient] = TVfunctional(Signal,Distribution,L,Kernel,RegularizationParameter)

Functional = 1/2*norm(Signal - Kernel*Distribution)^2 + RegularizationParameter^2*sum(sqrt((L*Distribution).^2 + 1e-24));
if nargout>1
    TVGradient = (L)'*diag(1./(sqrt((L*Distribution).^2 + 1e-24)))*(L)*Distribution;
    Gradient = Kernel'*(Kernel*Distribution - Signal) + RegularizationParameter^2*TVGradient;
end

end


function [Functional,Gradient] = TikhonovFunctional(Signal,Distribution,L,Kernel,RegularizationParameter)

Functional = 1/2*norm(Kernel*Distribution - Signal)^2 + RegularizationParameter^2*1/2*norm(L*Distribution)^2;
if nargout>1
    Gradient = Kernel'*(Kernel*Distribution - Signal) + RegularizationParameter^2*(L')*L*Distribution;
end

end

function [Functional,Gradient] = HuberFunctional(Signal,Distribution,L,Kernel,RegularizationParameter,HuberParameter)

HuberPenalty = (sqrt(1+(L*Distribution/HuberParameter).^2)-1);
HuberFunctional = sum(HuberPenalty);
%Get the Huber norm gradient
Functional = 1/2*norm(Signal - Kernel*Distribution)^2 + RegularizationParameter^2*HuberFunctional;
if nargout>1
    HuberGradient =  RegularizationParameter^2*(L/(HuberParameter))'*diag(1./sqrt((L*Distribution/HuberParameter).^2 + 1))*(L/(HuberParameter))*Distribution;
    Gradient = Kernel'*(Kernel*Distribution - Signal) + RegularizationParameter^2*HuberGradient;
end

end