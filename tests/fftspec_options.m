function [err,data,maxerr] = test(opt,olddata)

FreqAxis = linspace(-10,10,200);
t = linspace(0,1/mean(2*10)*100,100);

S = exp(-t).*cos(2*pi*5*t);
Spectrum = abs(fftshift(fft(S,2*length(S))));

[X1] = fftspec(t,S,'Type','abs','ZeroFilling',200);
[X2] = fftspec(t,S,'Type','real','Apodization',false);

err(1) = any(length(X1)~=200);
err(2) = any(length(X2)~=200);
err = any(err);
maxerr = length(X1 - 200);
data = [];

if opt.Display
  figure,clf
  subplot(121)
  plot(t,S)
  subplot(122)
  hold on
  plot(FreqAxis,Spectrum)
  plot(FrequencyAxis,Xabs)
end

end