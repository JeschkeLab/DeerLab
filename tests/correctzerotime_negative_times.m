function [pass,maxerr] = test(opt)

% Check that zero-time correction works with negative times

t = linspace(-2,5,400);
r = time2dist(t);
P = dd_onegaussian(r,[4,0.2]);
S = dipolarkernel(t,r)*P;
zt = abs(min(t));

[ct,czt] = correctzerotime(S,t+1);

% Pass 1: corrected time-axis is equal to original
pass(1) = all(abs(ct - t.') < 1e-10);
% Pass 2: zero-time is returned properly
pass(2) = abs(czt - zt) < 1e-10;

pass = all(pass);
 
maxerr = max(abs(ct - t.'));

if opt.Display
   plot(t,S,ct,S)
   xlabel('t [\mus]')
   ylabel('V(t)')
   grid on, axis tight, box on
   legend('Input','Output')
end

end