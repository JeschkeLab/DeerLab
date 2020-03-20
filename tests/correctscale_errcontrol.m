function [pass,maxerr] = test(opt)

% Check error control of correctscale() towards wrong inputs

t = linspace(-1,5,200);
r = time2dist(t);
P = dd_onegaussian(r,[5,0.2]);
B = bg_exp(t,0.3);
TrueOffset = 1e8;
V = dipolarsignal(t,r,P,'Moddepth',0.25,'Background',B,'Scale',TrueOffset,'Phase',pi/2);

%Pass 1: passing a complex valued signal
try
    correctscale(V,t);
    pass = false;
catch
    pass = true;
end
 
maxerr = NaN;


end