function [err,data] = test(opt,olddata)

currentpath = pwd;
cd(fileparts(mfilename('fullpath')))
cd ../functions/private


%======================================================
% Grain down function test
%======================================================

dt = 4;
Grain = dt:dt:300;
originalData = 2 - log10(Grain);
dtDown = 2;
GrainDown = min(Grain):dtDown:max(Grain);
grainDownData = 2 - log10(GrainDown);



[grainedOutput,grainedt] = changegrain(originalData,Grain',dtDown);


err(1) = any(abs(grainedOutput - grainDownData)>1e-2);
err(2) = any(abs(grainedt - GrainDown)>1e-10);

err = any(err);
data = [];

cd(currentpath)

end