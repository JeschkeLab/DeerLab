function [pass,maxerr] = test(opt)

rng(2)
t = linspace(-1,4,100);
S = dipolarsignal(t,3,'noiselevel',0.02);

level1 = noiselevel(S);
level2 = noiselevel(S);


err(1) = ~isequal(round(level1,5),round(level2,5));

pass = all(err);

maxerr = max(abs(level1 - level2));
 

end