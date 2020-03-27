function [pass,maxerr] = test(opt)

% Check indifference of fftspec() towards input dimensionality

t = linspace(0,5,100);
S = rand(1,100);

[nu1,spec1] = fftspec(t,S);
[nu2,spec2] = fftspec(t.',S);
[nu3,spec3] = fftspec(t,S.');
[nu4,spec4] = fftspec(t.',S.');

% Pass 1: all frequency axes are equal
pass(1) = isequal(nu1,nu2,nu3,nu4);
% Pass 2: all spectra are equal
pass(2) = isequal(spec1,spec2,spec3,spec4);
% Pass 1: all frequency axes are column vectors
pass(3) = iscolumn(nu1) & iscolumn(nu2) & iscolumn(nu3) & iscolumn(nu4);
% Pass 1: all spectra are column vectors
pass(4) = iscolumn(spec1) & iscolumn(spec2) & iscolumn(spec3) & iscolumn(spec4);

pass = all(pass);

maxerr = NaN;
 

end