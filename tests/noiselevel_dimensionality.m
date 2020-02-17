function [err,data,maxerr] = test(opt,olddata)

rng(2)
t = linspace(-1,4,100);
S = dipolarsignal(t,3,'noiselevel',0.02);

level1 = noiselevel(S);
level2 = noiselevel(S);


err(1) = ~isequal(round(level1,5),round(level2,5));

err = any(err);

maxerr = max(abs(level1 - level2));
data = [];

end