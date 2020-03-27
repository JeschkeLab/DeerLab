function [pass,maxerr] = test(opt)

% Check that overtones are properly generated in dipolar kernel with all methods

t = linspace(0,3,50); % us
r = time2dist(t); % nm
Tmix = 50; % us
T1 = 88; % us

coefficients = [0.50 0.29 0.20];

Kfresnel = dipolarkernel(t,r,'Method','fresnel','OvertoneCoeffs',coefficients);
Kgrid = dipolarkernel(t,r,'Method','grid','OvertoneCoeffs',coefficients);
Kint = dipolarkernel(t,r,'Method','grid','OvertoneCoeffs',coefficients);

delta1 = abs(Kfresnel - Kgrid);
delta2 = abs(Kfresnel - Kint);
delta3 = abs(Kint - Kgrid);

% Pass 1: fresnel and grid kernels are equal
pass(1) = all(delta1(:) < 1e-4);
% Pass 2: fresnel and integral kernels are equal
pass(2) = all(delta2(:) < 1e-4);
% Pass 3: integral and grid kernels are equal
pass(3) = all(delta3(:) < 1e-4);

pass = all(pass);

maxerr = max(delta1(:));

 

end