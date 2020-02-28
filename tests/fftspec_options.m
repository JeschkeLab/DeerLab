function [pass,maxerr] = test(opt)

%Check that fftspec() runs with all its options

nu = linspace(-10,10,200);
t = linspace(0,1/mean(2*10)*100,100);
S = exp(-t).*cos(2*pi*5*t);
[nu1,spec1] = fftspec(t,S,'Type','abs','ZeroFilling',400);
[nu2,spec2] = fftspec(t,S,'Type','real','Apodization',false);


% Pass 1: spectra has right size
pass(1) = length(spec1) == 400;
% Pass 2: frequency axis has right size
pass(2) = length(nu1) == 400;
% Pass 2: spectra has right size
pass(3) = length(spec2) == 200;
% Pass 2: frequency axis has right size
pass(4) = length(nu2) == 200;

pass = all(pass);

maxerr = NaN;
 

if opt.Display
    plot(nu1,spec1,nu2,spec2)
    xlabel('\nu [MHz]')
    ylabel('Intensity [a.u.]')
    grid on, axis tight, box on
end

end