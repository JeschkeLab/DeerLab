function [err,data] = test(opt,olddata)

%======================================================
% Zero-time correction function
%======================================================

originalTimeAxis = -40:40;
originalData = 1000 - (originalTimeAxis.^2);
TimeAxis = originalTimeAxis + abs(min(originalTimeAxis));
inputZeroTime = TimeAxis(41);

[correctedTimeAxis,outputZeroTime] = correctZeroTime(originalData,TimeAxis,inputZeroTime);


err(1) = any(abs(correctedTimeAxis - originalTimeAxis)>1e-10);
err(2) = abs(outputZeroTime - inputZeroTime)>1e-10;

err = any(err);
data = [];

end