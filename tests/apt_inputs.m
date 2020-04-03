function [pass,maxerr] = test(opt)

% Check that apt() works with defined input schemes

%Parameters
t = linspace(-0.6,2,200);
r = time2dist(t);
P = dd_gauss(r,[3,0.5]);

K = dipolarkernel(t,r);
S = K*P;

aptK = aptkernel(t);

%Pass default value or let function use it
[aptP1] = apt(S,aptK,0.05);
[aptP2,r] = apt(S,aptK);

%Pass if the signal is equal down to numerical error
pass = isequal(aptP1,aptP2);
 
maxerr = NaN;

%Plot if requested
if opt.Display 
    plot(r,aptP1,r,aptP2);
    xlabel('r [nm]')
    ylabel('P(r) [nm^{-1}]')
    legend('Input scheme 1','Input scheme 2')
    axis tight, grid on, box on
end


end