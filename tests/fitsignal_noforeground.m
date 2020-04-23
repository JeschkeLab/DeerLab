function [pass,maxerr] = test(opt)

% Check that fitsignal() can fit without any foreground

rng(1)

t = linspace(0,5,100);
r = linspace(2,6,200);
B = bg_exp(t,0.2);

[~,~,Bfit] = fitsignal(B,t,r,'none',@bg_exp,@exp_4pdeer);

% Pass 1: background is well fitted
pass = all(abs(B - Bfit) < 1e-8);

maxerr = max(abs(B - Bfit));

if opt.Display
   plot(t,B,'k.',t,Bfit,'r')
   legend('truth','fit')
   xlabel('t [\mus]')
   ylabel('B(t)')
   grid on, axis tight, box on
end

end