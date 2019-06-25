% Based on regularization tools by Per Christian Hansen
% see: http://www2.compute.dtu.dk/~pcha/Regutools/index.html
% as well as: http://ch.mathworks.com/matlabcentral/fileexchange/52-regtools

noiselev = 0.05; % noise level for tests, 0.001 to 0.1 is more or less fine

load pake_base_tikh_512 % kernel data with singular value decomposition
                    % see make_bas_Tikh for kernel generation
par = [2.5,0.4,6,0.75,0.5]; % Two Gaussians at 3 and 4.5 nm with widths of 0.2 and 0.4 nm same integral intensity
distr = Two_Gaussians(r,par); % make test tata

figure(1); clf; % plot test data distribution
plot(r,distr);
title('Input distance distribution');

ff = get_ff(r,distr,t,kernel,r,t); % simulate form factor

ff = ff+noiselev*randn(1,length(ff)); % add pseudo-noise

data = [1000*t' ff' zeros(size(ff'))]; % save test data 
save test.dat data -ascii

% Compute L curve, it's corner, and make regularization at the L curve
% corner and compute form factor ff2 corresponding to the Tikhonov solution
tic,
[corner_index,idx_AIC,idx_GCV,rho,eta,reg_param] = l_curve_mod(kernel,L,ff');
distr2 = tikhonov(kernel,L,ff(:),reg_param(corner_index));
rho1 = norm(kernel*distr2-ff(:));
eta1 = norm(L*distr2);
ff2 = get_ff(r,distr2,t,kernel,r,t);
toc

distr2p = distr2;
distr2p(distr2<0)=0;
alpha = reg_param(corner_index);
tic,
C = [kernel';alpha*L];
d = [ff';L*distr2p];
distr3 = lsqnonneg(C,d);
ff3 = get_ff(r,distr3,t,kernel,r,t);
toc,

% plot distance distribution from Tikhonov regularization
figure(2); clf; 
plot(r,distr2);
hold on;
plot(r,distr3,'k');
title('Tikhonov regularization result');

% plot input form factor and its fit
figure(3); clf;
plot(t,ff);
hold on;
plot(t,ff2,'r');
plot(t,ff3,'b');
title('Input form factor and fit');

% plot L curve and display corner localization
figure(4); clf;
plot(rho,eta,'k.');
hold on;
plot(rho(corner_index),eta(corner_index),'ro');
title(sprintf('L curve and corner localization at %12.1f',reg_param(corner_index)));

% plot difference between form factors obtained with unsigned and
% non-negative Tikhonov regularization
figure(5); clf;
plot(t,ff3-ff,'k');
hold on;
plot(t,ff3-ff2,'r');
title('Fit residual and difference between unconstrained reg. and reg. with P>0');