function [err,data,maxerr] = test(opt,olddata)

t = linspace(-1,5,200);
r= time2dist(t);
P = rd_onegaussian(r,[5,0.2]);
B = td_exp(t,0.3);
TrueOffset = 1e8;
V = dipolarsignal(t,r,P,'Moddepth',0.25,'Background',B,'Scale',TrueOffset);

[Vc,Offset] = correctscale(V,t);

err  = abs(Offset - TrueOffset)>1e-5;
data = [];
maxerr = max(abs(Offset - TrueOffset));

if opt.Display
    figure(8)
    plot(t,Vc)   
end

end