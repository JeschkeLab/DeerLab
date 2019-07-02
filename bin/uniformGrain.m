function output=uniformGrain(nonUniformGrain,data,uniformGrain),
%
% Computes a new vector defined on a non-uniform grain from the input 
% vector data to a new vector defined on a uniform grain defined by the
% input. 
%
% (Adapted from G. Jeschke by Luis Fabregas
%

% Standard resolution 2000 points


Grain = uniformGrain(2)-uniformGrain(1);
minGrain = min(nonUniformGrain);
if minGrain < min(uniformGrain)
    minGrain = min(uniformGrain); 
end
maxGrain = max(nonUniformGrain);
if maxGrain > max(uniformGrain)
    maxGrain=max(uniformGrain); 
end

minPoint = round(minGrain/Grain);
maxPoint = round(maxGrain/Grain);
rax = linspace(minPoint*Grain, maxPoint*Grain, maxPoint - minPoint + 1);
InterpData = interp1(nonUniformGrain, data, rax, 'pchip', 0);
output = 0*uniformGrain;
rbas = round(min(uniformGrain)/Grain);
output(minPoint - rbas + 1:maxPoint - rbas + 1) = InterpData;



