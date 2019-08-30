function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check kernel is constructed properly
%======================================================

Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,Dimension*TimeStep,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = rd_onegaussian(DistanceAxis,[3,0.5]);
Distribution = Distribution/sum(Distribution)/mean(diff(DistanceAxis));
Background = exp(-0.5*TimeAxis);
Kernel = dipolarkernel(TimeAxis,DistanceAxis);

Trace  = Kernel*Distribution;
Trace = (Trace + 2).*Background';
Background = Background*(1-1/Trace(1));
Trace = Trace/Trace(1);

KernelB = dipolarkernel(TimeAxis,DistanceAxis,Background,'KernelBType','full');
TraceB  = KernelB*Distribution;

err = any(abs(TraceB - Trace)>1e-10);
maxerr = max(abs(TraceB - Trace));
data = [];

if opt.Display
   figure(3)
   plot(TimeAxis,TraceB,'r',TimeAxis,Background,'r--',TimeAxis,Trace,'b')
   legend('truth','B','K*P')
end

end