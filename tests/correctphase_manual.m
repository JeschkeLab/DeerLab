function [pass,maxerr] = test(opt)

% Test that manual phase correction works

x = (1:100).';
phiIn = pi/4;
xphased = x.*exp(1i*phiIn);
[xRe,xIm,phiOut] = correctphase(xphased,phiIn);

% Pass 1: real part is equal to unphased input
pass(1) = any(abs(xRe - real(x)) < 1e-10);
% Pass 2: imaginary part is equal to unphased input
pass(2) = any(abs(xIm - imag(x)) < 1e-10);
% Pass 3: phase is returned properly
pass(3) = isequal(phiIn,phiOut);

pass = all(pass);

maxerr = max(abs(xIm - imag(x)));
 

end