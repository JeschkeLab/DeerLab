function [err,data] = test(opt,olddata)

%======================================================
% Grain up function test
%======================================================

TimeStep = 2;
Grain = TimeStep:TimeStep:100;
originalData = 2 - log10(Grain);
TimeStepUp = 8;
GrainUp = TimeStepUp:TimeStepUp:100;
grainUpData = 2 - log10(GrainUp);



[grainedOutput,grainedTimeAxis] = changegrain(originalData,Grain',TimeStepUp);


err(1) = any(abs(grainedOutput - grainUpData)>1e-10);
err(2) = any(abs(grainedTimeAxis - GrainUp)>1e-10);

err = any(err);
data = [];

end