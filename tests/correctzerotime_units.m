function [err,data,maxerr] = test(opt,olddata)

originalt = -5:0.5:80;
originalData = 1000 - (originalt.^2);
originalData = originalData/max(originalData);
t = originalt + abs(min(originalt));
inputZeroTime = abs(min(originalt));
%ns
correctedt1 = correctzerotime(originalData,t);
%us
t = t/1000;
correctedt2 = correctzerotime(originalData,t);

err = any(abs(correctedt1 - correctedt2*1000)>1e-7);
maxerr = max(abs(correctedt1 - correctedt2));
data = [];

end