function [pass,maxerr] = test(opt)

% Check that selectmodel() works with both LSQ and Chi^2 cost functionals

t = linspace(0,3,200);
r = linspace(2,6,100);
parIn = [3 0.3 5 0.3 0.5];
P = dd_gauss2(r,parIn);
K = dipolarkernel(t,r);
S = K*P;
models = {@dd_gauss,@dd_gauss2,@dd_gauss3};

%nm
[opt1,metric1] = selectmodel(models,S,r,K,'aicc','CostModel','chi2red');
[opt2,metric2] = selectmodel(models,S,r,K,'aicc','CostModel','lsq');

% Pass: both cost functionals lead to the same result
pass = opt1==opt2;
 
maxerr = NaN;

if opt.Display
    subplot(121)
    plot(1:3,metric1)
    set(gca,'xtick',[1 2 3],'xticklabel',{'One Gaussian','Two Gaussian','Three Gaussian'})
    xtickangle(45)
    ylabel('AICc (LSQ)')
    grid on, axis tight, box on
    
    subplot(122)
    plot(1:3,metric2)
    set(gca,'xtick',[1 2 3],'xticklabel',{'One Gaussian','Two Gaussian','Three Gaussian'})
    xtickangle(45)
    ylabel('AICc (\chi^2)')
    grid on, axis tight, box on
end

end

