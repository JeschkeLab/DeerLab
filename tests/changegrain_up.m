function [err,data] = test(opt,olddata)

currentpath = pwd;
cd(fileparts(mfilename('fullpath')))
cd ../functions/private


%======================================================
% Grain up function test
%======================================================

dt = 2;
Grain = dt:dt:100;
originalData = 2 - log10(Grain);
dtUp = 8;
GrainUp = dtUp:dtUp:100;
grainUpData = 2 - log10(GrainUp);



[grainedOutput,grainedt] = changegrain(originalData,Grain',dtUp);


err(1) = any(abs(grainedOutput - grainUpData)>1e-10);
err(2) = any(abs(grainedt - GrainUp)>1e-10);

err = any(err);
data = [];

cd(currentpath)


end