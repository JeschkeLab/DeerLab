function [pass,maxerr] = test(opt)

% Test correctphase() for 2D datasets with given phases

x = repmat((1:100).',1,20);
phases = linspace(pi/8,pi/2,20);
xphased = x.*exp(-1i*phases);

[xcorr,~,phasesfit] = correctphase(xphased,phases);

% Pass 1: the imagninary part if properly minimized
pass(1) = all(all(abs(imag(xcorr) - imag(x)) < 1e-10));
% Pass 2: the input phase is properly fitted
pass(2) = all(abs(phases - phasesfit) < 1e-10);
% Pass 3: output dimensions are maintained
pass(3) = isequal(size(xphased),size(xcorr));


pass = all(pass);

maxerr = max(max(abs(imag(xcorr) - imag(x))));
 

end