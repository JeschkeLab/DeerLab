function [pass,maxerr] = test(opt)

% Basic functionality test on dipolarbackground()

t = linspace(0,5,150);
kappa = 0.3;
lam = 0.5;

%Reference
B0 = bg_exp(t,lam*kappa);

%Output
Bmodel = @(t) bg_exp(t,kappa);
path(1,:) = [1-lam NaN];
path(2,:) = [lam 0];

B = dipolarbackground(t,path,Bmodel);

% Pass 1: background produced as expected
pass(1) = all(abs(B-B0) < 1e-8);
% Pass 2: background has expected maxima
[~,idx] = max(B);
pass(2) =  t(idx)==0;

pass = all(pass);

maxerr = max(abs(B-B0));

if opt.Display
   plot(t,B0,'k',t,B,'r')
   axis tight, grid on, box on
   xlabel('t [\mus]')
   ylabel('B(t)')
   legend('ref','out')
end

end
