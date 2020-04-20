function [pass,maxerr] = test(opt)

% Test dipolarbackground() returns correct multipathway background without
% any renormalization

t = linspace(0,5,150);
kappa = 0.3;
T0 = 2.5;
lams = [0.4 0.6 0.3];
lams = lams/sum(lams);

%Reference
B0 = bg_exp(t,lams(2)*kappa);
B0 = B0.*bg_exp(t-T0,lams(3)*kappa);

%Output
Bmodel = @(t) bg_exp(t,kappa);
path(1,:) = [lams(1) NaN];
path(2,:) = [lams(2) 0];
path(3,:) = [lams(3) T0];

B = dipolarbackground(t,path,Bmodel,'renormalize',false);

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
   ylim([0 1])
end

end
