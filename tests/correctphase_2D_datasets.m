function [pass,maxerr] = test(opt)

% Test correctphase() for 2D datasets

x = repmat((1:100).',1,20);
phases = mod(linspace(-3*pi/4,pi/2,20),pi);
xphased = x.*exp(1i*phases);

[xcorr,~,phasesfit] = correctphase(xphased);

% Pass 1: the imagninary part if properly minimized
pass(1) = all(all(abs(imag(xcorr) - imag(x)) < 1e-10));
% Pass 2: the input phase is properly fitted
pass(2) = all(abs(phases - phasesfit) < 1e-4);
% Pass 3: output dimensions are maintained
pass(3) = isequal(size(xphased),size(xcorr));

pass = all(pass);

maxerr = max(abs(phases + phasesfit));
 

end