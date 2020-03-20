function [pass,maxerr] = test(opt)

% Test that correctscale find accurately the scale in the presence of noise

t = linspace(-1,5,200);
r = time2dist(t);
P = dd_onegauss(r,[5,0.2]);
B = bg_exp(t,0.3);
scale = 1e8;
noiselevel = 0.05;
V = dipolarsignal(t,r,P,'Moddepth',0.25,'Background',B,...
    'Scale',scale,'noiselevel',noiselevel);

[Vc,scalefit] = correctscale(V,t);

% Pass: scale found accurately within noise level
pass  = abs(scalefit - scale)/scale < noiselevel;

maxerr = max(abs(scalefit - scale)/scale);

if opt.Display
    plot(t,V/scale,t,Vc)  
    xlabel('t [\mus]')
    ylabel('V(t) ')
    legend('Truth','Corrected')
    axis tight, grid on, box on
end

end