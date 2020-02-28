function [pass,maxerr] = test(opt)

t = linspace(-1,4,100);
S = dipolarsignal(t,3);

S1 = longpass(t,S);
S2 = longpass(t.',S);
S3 = longpass(t,S.');
S4 = longpass(t.',S.');


err(1) = ~isequal(S1,S2,S3,S4);
err(2) = ~iscolumn(S1) | ~iscolumn(S2) | ~iscolumn(S3) | ~iscolumn(S4);

pass = all(err);

maxerr = max(abs(S1 - S2));
 

end