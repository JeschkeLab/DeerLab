function [pass,maxerr] = test(opt)

% Check that ffstpec() runs correctly and returns the spectrum

nuref = linspace(-10,10,200);
t = linspace(0,1/mean(2*10)*100,100);

S = exp(-t).*cos(2*pi*5*t);
specRef = abs(fftshift(fft(S,2*length(S)))).';
[nu,spec] = fftspec(t,S,'Type','abs','Apodization',false);

error = abs(specRef - spec);

% Pass: the returned spectrum is equal to the reference
pass = all(error < 1e-10);

maxerr = max(error);
 
if opt.Display
    plot(nuref,specRef,nu,spec)
    xlabel('\nu [MHz]')
    ylabel('Intensity [a.u.]')
    grid on, axis tight, box on
    legend('reference','output')
end

end