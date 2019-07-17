function [DipEvoFcn] = supressGhostDistances(FormFactor,ModulationDepth,NRadicals)

Scaling = 1/(NRadicals - 1);

FormFactor = FormFactor.^Scaling;

DipEvoFcn = (FormFactor - (1 - ModulationDepth))/ModulationDepth;
DipEvoFcn = DipEvoFcn/DipEvoFcn(1);

end





