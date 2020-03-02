function [pass,maxerr] = test(opt)

% Check apodization has an effect on spectrum

t = linspace(0,1/mean(2*10)*100,100);
S = exp(-t).*cos(2*pi*5*t);
[nu1,spec1] = fftspec(t,S,'Type','abs','Apodization',false);
[nu2,spec2] = fftspec(t,S,'Type','abs','Apodization',true);

% Pass: spectra are different due to apodization
pass = ~isequal(spec1,spec2);

maxerr = NaN;

if opt.Display
    plot(nu1,spec1,nu2,spec2)
    xlabel('\nu [MHz]')
    ylabel('Intensity [a.u.]')
    grid on, axis tight, box on
end

end