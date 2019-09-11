%
% UNIFORMGRAIN Conversion from non-uniform grain to uniform grain
%
%   uD = UNIFORMGRAIN(nug,nuD,ug)
%   Computes a new data vector (uD) defined on a non-uniform grain (nug) 
%   from the input vector data  (nuD) to a new vector defined on a uniform
%   grain defined by the input (ug). 
%

% This file is a part of DeerAnalysis. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.



function output = uniformgrain(nonUniformGrain,data,uniformgrain)
%
% Computes a new vector defined on a non-uniform grain from the input 
% vector data to a new vector defined on a uniform grain defined by the
% input. 
%
% (Adapted from G. Jeschke by Luis Fabregas
%

% Standard resolution 2000 points


Grain = uniformgrain(2)-uniformgrain(1);
minGrain = min(nonUniformGrain);
if minGrain < min(uniformgrain)
    minGrain = min(uniformgrain); 
end
maxGrain = max(nonUniformGrain);
if maxGrain > max(uniformgrain)
    maxGrain=max(uniformgrain); 
end

minPoint = round(minGrain/Grain);
maxPoint = round(maxGrain/Grain);
rax = linspace(minPoint*Grain, maxPoint*Grain, maxPoint - minPoint + 1);
InterpData = interp1(nonUniformGrain, data, rax, 'pchip', 0);
output = 0*uniformgrain;
rbas = round(min(uniformgrain)/Grain);
output(minPoint - rbas + 1:maxPoint - rbas + 1) = InterpData;
output = output/sum(output)/mean(diff(uniformgrain));



