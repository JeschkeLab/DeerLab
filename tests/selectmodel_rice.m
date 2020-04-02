function [pass,maxerr] = test(opt)

% Test if selectmethod can identify that the optimal method is a two
% Rician model as given as the input signal

t = linspace(0,3,200);
r = linspace(2,6,100);
parIn = [3 0.3 5 0.3 0.5];
P = dd_tworice(r,parIn);
K = dipolarkernel(t,r);
S = K*P;
models = {@dd_onerice,@dd_tworice,@dd_threerice};

[optimum,metric] = selectmodel(models,S,r,K,'aicc','Solver','lsqnonlin');

% Pass: the optimal model is found
pass = optimum==2;
 
maxerr = NaN;


if opt.Display
    plot(1:3,metric)
    set(gca,'xtick',[1 2 3],'xticklabel',{'One Rice','Two Rice','Three Rice'})
    xtickangle(45)
    ylabel('AICc')
    grid on, axis tight, box on
end

end

