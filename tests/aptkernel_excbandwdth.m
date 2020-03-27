function [pass,maxerr] = test(opt)

% Check that aptkernel() works with limited excitation bandwidth

t = linspace(0,4,150);
K1 = aptkernel(t,'ExcitationBandwidth',1e8);
K2 = aptkernel(t);

% Pass: all elements in kernel base are numerically equal
pass = all(all(abs(K1.Base-K2.Base) < 1e-10));

maxerr = max(max(abs(K1.Base-K2.Base)));
 

end