function [err,data,maxerr] = test(opt,olddata)

t = linspace(-1,5,200);
r= time2dist(t);
P = rd_onegaussian(r,[5,0.2]);
B = td_exp(t,0.3);
TrueOffset = 1e8;
V = dipolarsignal(t,r,P,'Moddepth',0.25,'Background',B,'Offset',TrueOffset);

%us
[Vc1] = correctscale(V,t);
%ns
t = t*1000;
[Vc2] = correctscale(V,t);

err  = abs(Vc1 - Vc2)>1e-10;
data = [];
maxerr = max(abs(Vc1 - Vc2));

if opt.Display
    figure(8)
    plot(t,Vc)
end

end