function [pass,maxerr] = test(opt)

% Compare excluded-volume with pure homogeneous background

c = 100; % uM
lambda = 0.6;
R = 0.1; % nm

t = linspace(-1,5,301);
 
B1 = bg_hom3d(t,c,lambda);
B2 = bg_hom3dex(t,[c R],lambda);

maxerr = max(abs(B1-B2));

pass = maxerr<1e-6;

end
