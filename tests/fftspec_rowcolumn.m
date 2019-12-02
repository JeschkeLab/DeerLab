function [err,data,maxerr] = test(opt,olddata)

FreqAxis = linspace(-10,10,200);
t = linspace(0,1/mean(2*10)*100,100);

S = exp(-t).*cos(2*pi*5*t);
Spectrum = abs(fftshift(fft(S,2*length(S))));
[nu,out1] = fftspec(t,S,'Type','abs','Apodization',false);
[nu,out2] = fftspec(t,S.','Type','abs','Apodization',false);

error = abs(out1 - out2);
err(1) = any(error>1e-10);

[nu,out1] = fftspec(t,S,'Type','abs','Apodization',false);
[nu,out2] = fftspec(t.',S,'Type','abs','Apodization',false);

error = abs(out1 - out2);
err(2) = any(error>1e-10);

err = any(err);

maxerr = max(error);
data = [];

end