function [err,data] = test(opt,olddata)

%======================================================
% Grain up function test
%======================================================

FreqAxis = 5:0.5:100;
nonUniformAxis = FreqAxis.^(1/3);
uniformAxis = linspace(min(nonUniformAxis),max(nonUniformAxis),length(nonUniformAxis));

nonUniformGaussian = gaussian(nonUniformAxis,3,0.25);
uniformGaussian = gaussian(uniformAxis,3,0.25);

[uniformGaussianOut] = uniformGrain(nonUniformAxis,nonUniformGaussian,uniformAxis);


err(1) = any(abs(uniformGaussianOut - uniformGaussian)>1e-1);

err = any(err);
data = [];

if opt.Display
   figure,hold on
   plot(uniformAxis,uniformGaussian) 
   plot(uniformAxis,uniformGaussianOut) 
end

end