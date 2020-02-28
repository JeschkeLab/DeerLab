function [pass,maxerr] = test(opt)

% Test that correctzerotime works with defaults


originalt = -5:0.5:80;
originalData = 1000 - (originalt.^2);
originalData = originalData/max(originalData);
t = originalt + abs(min(originalt));
inputZeroTime = abs(min(originalt));

[correctedt,outputZeroTime] = correctzerotime(originalData,t);

% Pass 1: corrected axis is equal to input
pass(1) = all(abs(correctedt - originalt.') < 1e-10);
% Pass 2: the zero time is accurately estimated
pass(2) = abs(outputZeroTime' - inputZeroTime) < 1e-10;

pass = all(pass);

maxerr = max(abs(correctedt - originalt.'));


end