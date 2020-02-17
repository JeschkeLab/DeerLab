function [err,data,maxerr] = test(opt,olddata)

t = linspace(0,5,100);
S = rand(1,100);

[nu1,spec1] = fftspec(t,S);
[nu2,spec2] = fftspec(t.',S);
[nu3,spec3] = fftspec(t,S.');
[nu4,spec4] = fftspec(t.',S.');


err(1) = ~isequal(nu1,nu2,nu3,nu4);
err(2) = ~isequal(spec1,spec2,spec3,spec4);
err(3) = ~iscolumn(nu1) | ~iscolumn(nu2) | ~iscolumn(nu3) | ~iscolumn(nu4);
err(4) = ~iscolumn(spec1) | ~iscolumn(spec2) | ~iscolumn(spec3) | ~iscolumn(spec4);

err = any(err);

maxerr = max(abs(spec1 - spec2));
data = [];

end