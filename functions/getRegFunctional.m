function functionHandle = getRegFunctional(Method,Signal,L,Kernel,RegularizationParameter)

switch Method
    case 'tv'
    functionHandle = @(Distribution)TVfunctional(Signal,Distribution,L,Kernel,RegularizationParameter);
    case 'tikhonov'
    functionHandle = @(Distribution)TikhonovFunctional(Signal,Distribution,L,Kernel,RegularizationParameter);
end

end

function [Functional,Gradient] = TVfunctional(Signal,Distribution,L,Kernel,RegularizationParameter)
Functional = 1/2*norm(Signal - Kernel*Distribution)^2 + RegularizationParameter^2*sum(sqrt((L*Distribution).^2 + 1e-24));
TVGradient = (L)'*diag(1./(sqrt((L*Distribution).^2 + 1e-24)))*(L)*Distribution;
Gradient = Kernel'*(Kernel*Distribution - Signal) + RegularizationParameter^2*TVGradient;
end


function [Functional,Gradient] = TikhonovFunctional(Signal,Distribution,L,Kernel,RegularizationParameter)
Functional = 1/2*norm(Kernel*Distribution - Signal)^2 + RegularizationParameter^2*1/2*norm(L*Distribution)^2;
Gradient = Kernel'*(Kernel*Distribution - Signal) + RegularizationParameter^2*(L')*L*Distribution;
end

