function [err,data] = test(opt,olddata)

%======================================================
% Grain down function test
%======================================================

TimeStep = 4;
Grain = TimeStep:TimeStep:300;
originalData = 2 - log10(Grain);
TimeStepDown = 2;
GrainDown = min(Grain):TimeStepDown:max(Grain);
grainDownData = 2 - log10(GrainDown);



[grainedOutput,grainedTimeAxis] = changeGrain(originalData,Grain',TimeStepDown);


err(1) = any(abs(grainedOutput - grainDownData)>1e-2);
err(2) = any(abs(grainedTimeAxis - GrainDown)>1e-10);

err = any(err);
data = [];

end