function [pass,maxerr] = test(opt)

% Test imaginary offset fitting of correctphase()

x = (1:100).';
phiIn = pi/4;
ImOffsetIn = 4i;
xphased = x.*exp(1i*phiIn) + ImOffsetIn;
fitImOffset = true;

[xcorr,~,phiOut,ImOffsetOut] = correctphase(xphased,[],fitImOffset);

% Pass 1: the imaginary part is properly minimized
pass(1) = all(abs(imag(xcorr) - imag(x)) < 1e-10);
% Pass 2: the input phase is properly fitted
pass(2) = abs(phiIn - phiOut) < 1e-4;
% Pass 3: the input imaginary offset is properly fitted
pass(3) = abs(ImOffsetIn - ImOffsetOut) < 1e-4;

pass = all(pass);

maxerr = max(abs(ImOffsetIn - ImOffsetOut));
 
end