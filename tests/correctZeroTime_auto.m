function [err,data] = test(opt,olddata)

%======================================================
% Zero-time correction function
%======================================================

originalTimeAxis = -5:0.5:80;
originalData = 1000 - (originalTimeAxis.^2);
originalData = originalData/max(originalData);
TimeAxis = originalTimeAxis + abs(min(originalTimeAxis));
inputZeroTime = abs(min(originalTimeAxis));

[correctedSignal,correctedTimeAxis,outputZeroTime] = correctzerotime(originalData,TimeAxis);


err(1) = any(abs(correctedTimeAxis - originalTimeAxis)>1e-10);
err(2) = any(abs(correctedSignal' - originalData)>1e-10);
err(3) = abs(outputZeroTime' - inputZeroTime)>1e-10;

err = any(err);
data = [];

end