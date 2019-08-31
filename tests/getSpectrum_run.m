function [err,data,maxerr] = test(opt,olddata)

FreqAxis = linspace(-10,10,200);
t = linspace(0,1/mean(2*10)*100,100);

S = exp(-t).*cos(2*pi*5*t);
Spectrum = abs(fftshift(fft(S,2*length(S))));
Output = fftspec(t,S);
[FrequencyAxis,~] = fftspec(t,S,'Type','abs');

error = abs(Spectrum - Output);
err = any(error>1e-10);
maxerr = max(error);
data = [];

if opt.Display
  figure,clf
  subplot(121)
  plot(t,S)
  subplot(122)
  hold on
  plot(FreqAxis,Spectrum)
  plot(FrequencyAxis,Output)
end

end