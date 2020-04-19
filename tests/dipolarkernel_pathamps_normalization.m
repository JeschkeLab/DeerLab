function [pass,maxerr] = test(opt)

% Check that dipolarkernel() works with arbitrary pathway amplitudes
t = linspace(-5,5,200);
r = linspace(0,8,200);
pathway(1,:) = [0.2 NaN];
pathway(2,:) = [0.8 0];
pathway(3,:) = [0.5 3];

K = dipolarkernel(t,r,pathway)/mean(diff(r));
Strace = K(:,100);

% Pass: there are no crashes and the kernel is properly normalized
pass = round(max(Strace),2)==1;

maxerr = NaN;

end