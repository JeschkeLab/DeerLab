function [err,data,maxerr] = test(opt,olddata)

S = rand(1,100) + 1i*rand(1,100);

[Vr1,Vi1] = correctphase(S);
[Vr2,Vi2] = correctphase(S.');


err(1) = ~isequal(Vr1,Vr2);
err(2) = ~isequal(Vi1,Vi2);
err(3) = ~iscolumn(Vr1) | ~iscolumn(Vr1);
err(4) = ~iscolumn(Vi1) | ~iscolumn(Vi1);

err = any(err);
maxerr = max(abs(Vr1 - Vr2));
data = [];

end