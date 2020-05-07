function [pass,maxerr] = test(opt)

% Check error control of correctphase() towards wrong inputs

V0 = 1:100;
inputPhase = pi/4;
V = V0.*exp(-1i*inputPhase);

%Pass 1: no inputs arguments
try
    correctphase;
    pass(1) = false;
catch
    pass(1) = true;
end

% Pass 2: wrong input arguments
try
    correctphase(V,inputPhase,false,'wrongargument');
    pass(2) = false;
catch
    pass(2) = true;
end

pass = all(pass);

maxerr = NaN;

end