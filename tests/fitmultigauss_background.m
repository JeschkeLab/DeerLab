function [pass,maxerr] = test(opt)

% Check that fitmultigauss() can fit the background and modulation depth

t = linspace(0,5,100);
r = time2dist(t);
InputParam = [4 0.2 4 1 3 0.4 0.4 0.4];
P = dd_threegaussian(r,InputParam);
B = bg_exp(t,0.55);
V = dipolarsignal(t,r,P,'moddepth',0.5,'background',B);
FitP = fitmultigauss(V,t,r,5,'aicc','background',@bg_exp,'Upper',[6 1 0.9 1],'Lower',[1 0.1 0.2 0.01]);

% Pass: distribution is well fitted
pass = all(abs(FitP - P) < 8e-1);

maxerr = max(abs(FitP - P));

if opt.Display
   plot(r,P,'k',r,FitP,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end