function [pass,maxerr] = test(opt)

% Check indifference of noiselevel() towards input dimensionality

rng(1)
t = linspace(-1,4,100);
S = dipolarsignal(t,3,'noiselevel',0.02);

level1 = noiselevel(S);
level2 = noiselevel(S.');

% Pass: noiselevels are equal
pass = isequal(round(level1,5),round(level2,5));

maxerr = NaN;
 

end