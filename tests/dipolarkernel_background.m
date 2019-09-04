function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check kernel is constructed properly
%======================================================

Dimension = 200;
dt = 0.008;
t = linspace(0,Dimension*dt,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);
P = P/sum(P)/mean(diff(r));
B = exp(-0.5*t);
K = dipolarkernel(t,r);

Trace  = K*P;
Trace = (Trace + 2).*B';
ModDepth = 1/Trace(1);
Trace = Trace/Trace(1);

KB = dipolarkernel(t,r,ModDepth,B);
TraceB  = KB*P;

err = any(abs(TraceB - Trace)>1e-10);
maxerr = max(abs(TraceB - Trace));
data = [];

if opt.Display
   figure(3)
   plot(t,TraceB,'r',t,B,'r--',t,Trace,'b')
   legend('truth','B','K*P')
end

end