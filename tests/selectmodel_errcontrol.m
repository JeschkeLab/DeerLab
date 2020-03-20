function [pass,maxerr] = test(opt)

% Check error control of regoperator() towards wrong inputs

t = linspace(0,5,80);
S = dipolarsignal(t,3);
r = linspace(1,6,50);
K = dipolarkernel(t,r);
Models = {@dd_onegauss,@dd_twogauss,@dd_threegauss};

% Pass 1: invalid selection method
try
    selectmodel(Models,S,r,K,'gsad'); 
    pass(1) = false;
catch
    pass(1) = true;
end

% Pass 2: not enough input arguments
try
    selectmodel(Models); 
    pass(2) = false;
catch
    pass(2) = true;
end

% Pass 3: not enough upper boundaries
try
    selectmodel(Models,S,r,K,'aic','upper',[1 2 3]); 
    pass(3) = false;
catch
    pass(3) = true;
end

% Pass 4: not enough lower boundaries
try
    selectmodel(Models,S,r,K,'aic','lower',[1 2 3]); 
    pass(4) = false;
catch
    pass(4) = true;
end

% Pass 5: too many lower boundaries
try
    selectmodel(Models,K,'aic','lower',{[1 2 3],[1 2 3],[1 3 4],[1 4 5]}); 
    pass(5) = false;
catch
    pass(5) = true;
end

% Pass 6: option string not defined
try
    selectmodel(Models,S,r,K,'aic',[1 2 3]); 
    pass(6) = false;
catch
    pass(6) = true;
end

pass = all(pass);

maxerr = NaN;
 

end

