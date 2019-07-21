function [err,data,maxerr] = test(opt,olddata)

FreqAxis = linspace(-10,10,200);
TimeAxis = linspace(0,1/mean(2*10)*100,100);

Signal = exp(-TimeAxis).*cos(2*pi*5*TimeAxis);
Spectrum = abs(fftshift(fft(Signal,2*length(Signal))));
Output = getSpectrum(TimeAxis,Signal);
[FrequencyAxis,~] = getSpectrum(TimeAxis,Signal,'Type','abs');

error = abs(Spectrum - Output);
err = any(error>1e-10);
maxerr = max(error);
data = [];

if opt.Display
  figure,clf
  subplot(121)
  plot(TimeAxis,Signal)
  subplot(122)
  hold on
  plot(FreqAxis,Spectrum)
  plot(FrequencyAxis,Output)
end

end