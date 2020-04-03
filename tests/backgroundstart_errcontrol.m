function [pass,maxerr] = test(opt)

% Check error control of backgroundstart() towards wrong inputs

t = linspace(0,4,150);
r = time2dist(t);
P = dd_gauss(r,[3,0.5]);
B = bg_exp(t,0.5);
lam0 = 0.5;
S = dipolarsignal(t,r,P,'background',B,'moddepth',lam0);

% Pass 1: forgetting the background model
try
    backgroundstart(S,t);
    pass(1) = false;
catch
    pass(1) = true;
end

% Pass 2: passing a complex signal
S2 = dipolarsignal(t,r,P,'background',B,'moddepth',lam0,'phase',pi/2);
try
    backgroundstart(S2,t,@bg_poly1);
    pass(2) = false;
catch
    pass(2) = true;
end

% Pass 3: passing a non function-handle as model
model = rand(2,2);
try
    backgroundstart(S,t,model);
    pass(3) = false;
catch
    pass(3) = true;
end

% Pass 3: passing relative start/end positions in inverse order
try
    backgroundstart(S,t,@bg_poly1,'RelSearchStart',0.9,'RelSearchEnd',0.1);
    pass(4) = false;
catch
    pass(4) = true;
end

pass = all(pass);
 
maxerr = NaN;

end