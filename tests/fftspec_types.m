function [pass,maxerr] = test(opt)

% Check that the different spectra types are computed correctly 

t = linspace(0,1/mean(2*10)*100,100);
S = exp(-t).*cos(2*pi*5*t);

[nu,specAbs] = fftspec(t,S,'Type','abs');
specRe = fftspec(t,S,'Type','real');
specIm = fftspec(t,S,'Type','imag');

error = abs(sqrt(specRe.^2 + specIm.^2) - specAbs);

% Pass: absolute spectrum is can be computed from the real/imag spectra
pass = all(error < 1e-10);

maxerr = max(error);
 
if opt.Display
    plot(nu,specAbs,nu,specRe,nu,specIm)
    xlabel('\nu [MHz]')
    ylabel('Intensity [a.u.]')
    grid on, axis tight, box on
    legend('reference','output')
end

end