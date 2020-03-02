function [pass,maxerr] = test(opt)

% Test that correctscale find accurately the scale in abscence of noise

t = linspace(-1,5,200);
r = time2dist(t);
P = rd_onegaussian(r,[5,0.2]);
B = td_exp(t,0.3);
scale = 1e8;
V = dipolarsignal(t,r,P,'Moddepth',0.25,'Background',B,'Scale',scale);

[Vc,scalefit] = correctscale(V,t);

% Pass: scale found accurately
pass  = abs(scalefit - scale) < 1e-7;
 
maxerr = max(abs(scalefit - scale));

if opt.Display
    plot(t,V/scale,t,Vc)  
    xlabel('t [\mus]')
    ylabel('V(t) ')
    legend('Truth','Corrected')
    axis tight, grid on, box on
end

end