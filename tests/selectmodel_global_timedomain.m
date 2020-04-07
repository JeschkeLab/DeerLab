function [pass,maxerr] = test(opt)

% Test if selectmethod can identify the optimal background model


kappa = 0.2;
d = 2.5;

t1 = linspace(0,3,200);
B1 = bg_strexp(t1,[kappa d]);

t2 = linspace(0,3,200);
B2 = bg_strexp(t2,[kappa d]);

t3 = linspace(0,3,200);
B3 = bg_strexp(t3,[kappa d]);

models = {@bg_exp,@bg_strexp,@bg_prodstrexp};

[optimum,metric] = selectmodel(models,{B1,B2,B3},{t1,t2,t3},'aicc');

% Pass: all selection methods find the optimal model
pass = optimum==2;
 
maxerr = NaN;

if opt.Display
    plot(1:3,metric)
    set(gca,'xtick',[1 2 3],'xticklabel',{'Exponential','Stretched Exponential','Product of stretched exponentials'})
    xtickangle(45)
    ylabel('AICc')
    grid on, axis tight, box on
end

end

