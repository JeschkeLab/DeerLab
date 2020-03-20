function [pass,maxerr] = test(opt)

% Test if selectmethod() can identify that the optimal method is a two
% Gaussian model as given as the input signal

t = linspace(0,3,200);
r = linspace(2,6,100);
parIn = [3 0.3 5 0.3 0.5];
P = dd_twogaussian(r,parIn);
K = dipolarkernel(t,r);
S = K*P;
models = {@dd_onegaussian,@dd_twogaussian,@dd_threegaussian};

[optimum1,metrics1] = selectmodel(models,S,r,K,{'aic','aicc','bic','rmsd'});
[optimum2,metrics2] = selectmodel(models,S,r,K,'all');

% Pass 1: both input schemes return the same optimums
pass(1) = isequal(optimum1,optimum2);
% Pass 2: all selection functionals find the same optimal model
pass(2) = isequal(optimum1(1),optimum1(2),optimum1(3));
% Pass 3: both input schemes return the same functionals
pass(3) = isequal(metrics1,metrics2);

pass = all(pass);
 
maxerr = NaN;


if opt.Display
    subplot(221)
    plot(1:3,metrics1{1})
    set(gca,'xtick',[1 2 3],'xticklabel',{'One Gaussian','Two Gaussian','Three Gaussian'})
    xtickangle(45)
    ylabel('AIC')
    grid on, axis tight, box on
    
    subplot(222)
    plot(1:3,metrics1{2})
    set(gca,'xtick',[1 2 3],'xticklabel',{'One Gaussian','Two Gaussian','Three Gaussian'})
    xtickangle(45)
    ylabel('AICc')
    grid on, axis tight, box on
    
    subplot(223)
    plot(1:3,metrics1{3})
    set(gca,'xtick',[1 2 3],'xticklabel',{'One Gaussian','Two Gaussian','Three Gaussian'})
    xtickangle(45)
    ylabel('BIC')
    grid on, axis tight, box on
    
    subplot(224)
    plot(1:3,metrics1{4})
    set(gca,'xtick',[1 2 3],'xticklabel',{'One Gaussian','Two Gaussian','Three Gaussian'})
    xtickangle(45)
    ylabel('RMSD')
    grid on, axis tight, box on
end

end

