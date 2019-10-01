function [err,data,maxerr] = test(opt,olddata)

t = linspace(-1,5,200);
r= time2dist(t);
P = rd_onegaussian(r,[5,0.2]);
B = td_exp(t,0.3);
TrueOffset = 1e8;
V = dipolarsignal(t,r,P,'Moddepth',0.25,'Background',B,'Offset',TrueOffset,'Phase',pi/2);

try
    Vc = correctscale(V,t);
    err = true;
catch
    err = false;
end
data = [];
maxerr = 0;


end