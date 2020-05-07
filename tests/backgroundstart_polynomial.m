function [pass,maxerr] = test(opt)

% Check that backgroundstart() works with polynomial background signals

t = linspace(0,3.2,200);
r = time2dist(t);
S = dipolarkernel(t,r)*dd_gauss(r,[3,0.5]);
lam = 0.5;
F = (1 - lam) + lam*S;
B = 1  - lam*t.';
V = F.*B;
t0 = backgroundstart(V,t,@bg_poly1);
[Bfit,lam] = fitbackground(V,t,@bg_poly1,t0);
Bfit = Bfit*(1-lam);
B = B*(1-lam);

%Pass: the background fitted using the optimized start fits well
pass = abs(B - Bfit) < 7e-4;
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