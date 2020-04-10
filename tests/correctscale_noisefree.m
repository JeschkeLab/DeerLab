function [pass,maxerr] = test(opt)

% Test that correctscale find accurately the scale in abscence of noise

t = linspace(-1,5,300);
r = time2dist(t);
P = dd_gauss(r,[5,0.2]);
B = bg_exp(t,0.3);
scale = 1e8;
V = dipolarsignal(t,r,P,0.25,B,'Scale',scale);

[Vc,scalefit] = correctscale(V,t);

% Pass: scale found accurately
pass  = abs(scalefit/scale - 1) < 1e-12;
 
maxerr = max(abs(scalefit - scale));

if opt.Display
    plot(t,V/scale,t,Vc)  
    xlabel('t [\mus]')
    ylabel('V(t) ')
    legend('Truth','Corrected')
    axis tight, grid on, box on
end

end