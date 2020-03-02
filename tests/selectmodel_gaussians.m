function [pass,maxerr] = test(opt)

% Test if selectmethod can identify that the optimal method is a two
% Gaussian model as given as the input signal


t = linspace(0,3,200);
r = linspace(2,6,100);
parIn = [3 0.3 5 0.3 0.5];
P = rd_twogaussian(r,parIn);
K = dipolarkernel(t,r);
S = K*P;
models = {@rd_onegaussian,@rd_twogaussian,@rd_threegaussian};

[optimum1,metric1] = selectmodel(models,S,r,K,'aicc');
[optimum2,metric2] = selectmodel(models,S,r,K,'aic');
[optimum3,metric3] = selectmodel(models,S,r,K,'bic');
[optimum4,metric4] = selectmodel(models,S,r,K,'rmsd');

% Pass: all selection methods find the optimal model
pass = optimum1==2 & optimum2==2 & optimum3==2 & optimum4==2;
 
maxerr = NaN;

if opt.Display
    subplot(221)
    plot(1:3,metric1)
    set(gca,'xtick',[1 2 3],'xticklabel',{'One Gaussian','Two Gaussian','Three Gaussian'})
    xtickangle(45)
    ylabel('AIC')
    grid on, axis tight, box on
    
    subplot(222)
    plot(1:3,metric2)
    set(gca,'xtick',[1 2 3],'xticklabel',{'One Gaussian','Two Gaussian','Three Gaussian'})
    xtickangle(45)
    ylabel('AICc')
    grid on, axis tight, box on
    
    subplot(223)
    plot(1:3,metric3)
    set(gca,'xtick',[1 2 3],'xticklabel',{'One Gaussian','Two Gaussian','Three Gaussian'})
    xtickangle(45)
    ylabel('BIC')
    grid on, axis tight, box on
    
    subplot(224)
    plot(1:3,metric4)
    set(gca,'xtick',[1 2 3],'xticklabel',{'One Gaussian','Two Gaussian','Three Gaussian'})
    xtickangle(45)
    ylabel('RMSD')
    grid on, axis tight, box on
end

end

