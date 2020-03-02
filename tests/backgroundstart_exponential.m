function [pass,maxerr] = test(opt)

% Check that backgroundstart() works with exponential background signals

t = linspace(0,3.2,200);
r = time2dist(t);
S = dipolarkernel(t,r)*rd_onegaussian(r,[3,0.5]);
B = td_exp(t,0.5);
lam0 = 0.5;
F = (1 - lam0) + lam0*S;
V = F.*B;

t0 = backgroundstart(V,t,@td_exp);
[Bfit] = fitbackground(V,t,@td_exp,t0);

%Pass: the background fitted using the optimized start fits well
pass = abs(B - Bfit) < 1e-4;
maxerr = max(abs(B - Bfit));
 
%Plot if requested
if opt.Display
    plot(t,B,t,Bfit)
    xlabel('t [\mus]')
    ylabel('B(t) ')
    legend('Truth','Fit')
    axis tight, grid on, box on
end

end