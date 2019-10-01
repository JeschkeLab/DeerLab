function [err,data,maxerr] = test(opt,olddata)

V0 = 1:100;
inputPhase = pi/4;

V = V0.*exp(-1i*inputPhase);

try
    correctphase;
    err = true;
catch
   err = false; 
end

try
    correctphase(V,inputPhase,false,'wrongargument');
    err = true;
catch
   err = false; 
end

data = [];
maxerr = 0;
end