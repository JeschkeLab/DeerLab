function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Zero-time correction function
%======================================================

originalt = -5:0.5:80;
originalData = 1000 - (originalt.^2);
originalData = originalData/max(originalData);
t = originalt + abs(min(originalt));
inputZeroTime = abs(min(originalt));

[correctedt,outputZeroTime] = correctzerotime(originalData,t);


err(1) = any(abs(correctedt - originalt)>1e-10);
err(2) = abs(outputZeroTime' - inputZeroTime)>1e-10;

maxerr = max(abs(correctedt - originalt));
err = any(err);
data = [];

end