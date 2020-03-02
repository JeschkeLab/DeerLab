function [pass,maxerr] = test(opt)

% Check that longpass() suppresses a single large frequency (short distance)

t = linspace(0,3,200);
r = linspace(0,6,200);
Ptrue = rd_onegaussian(r,[4,0.5]);
K = dipolarkernel(t,r);
S = K*Ptrue;

rart = 1.5;
freq = 52.04/(rart^3);
mod = 0.5*(2 + exp(-7*t).*cos(2*pi*freq*t));
S = S.*mod';
Sfilt = longpass(t,S,2);

Pfilt = fitregmodel(Sfilt,K,r,'tikhonov','aic','Solver','fnnls');
Praw = fitregmodel(Sfilt,K,r,'tikhonov','aic','Solver','fnnls');

error = abs(Pfilt - Ptrue);
%Pass: the large frequency is properly suppressed
pass = any(error < 2e-1);

maxerr = max(error);
 

if opt.Display
    
    subplot(131)
    plot(t,S,t,Sfilt)
    legend('raw','filtered')
    xlabel('t [\mus]')
    ylabel('V(t)')
    grid on, axis tight, box on
    
    subplot(132)
    plot(r,Ptrue,r,Praw,r,Pfilt)
    legend('truth','raw','filtered')
    xlabel('r [nm]')
    ylabel('P(r) [nm^{-1}]')
    grid on, axis tight, box on

    subplot(133),hold on
    [nu,specraw] = fftspec(t(t>=0),S(t>=0));
    specfilt = fftspec(t(t>=0),Sfilt(t>=0));
    plot(nu,specraw,nu,specfilt)
    legend('raw','filtered')
    xlabel('\nu [MHz]')
    ylabel('V(t)')
    grid on, axis tight, box on
    
end

end