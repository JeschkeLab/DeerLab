function [pass,maxerr] = test(opt)

% Check error control of aptkernel() towards wrong inputs

t = linspace(0,3,200);

% Pass 1: giving a string as excitation bandwidth
try
    K = aptkernel(t,'ExcitationBandwidth','a');
    pass = false;
catch
    pass = true;
end

maxerr = NaN;
 

end