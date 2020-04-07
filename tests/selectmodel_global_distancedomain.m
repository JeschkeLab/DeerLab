function [pass,maxerr] = test(opt)

% Test if selectmethod can identify that the optimal method is a two
% Gaussian model given multiple signals using global fitting


r = linspace(2,6,100);
parIn = [3 0.3 5 0.3 0.5];
P = dd_gauss2(r,parIn);

t1 = linspace(0,3,200);
K1 = dipolarkernel(t1,r);
S1 = K1*P;

t2 = linspace(0,1,200);
K2 = dipolarkernel(t2,r);
S2 = K2*P;

t3 = linspace(0,6,200);
K3 = dipolarkernel(t3,r);
S3 = K3*P;


models = {@dd_gauss,@dd_gauss2,@dd_gauss3};

[optimum,metric] = selectmodel(models,{S1,S2,S3},r,{K1,K2,K3},'aicc');

% Pass: all selection methods find the optimal model
pass = optimum==2;
 
maxerr = NaN;

if opt.Display
    plot(1:3,metric)
    set(gca,'xtick',[1 2 3],'xticklabel',{'One Gaussian','Two Gaussian','Three Gaussian'})
    xtickangle(45)
    ylabel('AICc')
    grid on, axis tight, box on
end

end

