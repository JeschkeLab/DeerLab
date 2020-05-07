function [pass,maxerr] = test(opt)

% Test if selectmethod can identify an optimal background model

t = linspace(0,3,200);
r = time2dist(t);
B = bg_strexp(t,[0.2 4]);
models = {@bg_exp,@bg_strexp,@bg_prodstrexp};

[optimum1,metric] = selectmodel(models,B,t,'aicc');
optimum2 = selectmodel(models,B,t,'aic');
optimum3 = selectmodel(models,B,t,'bic');

% Pass: the different functional find the same optimum
pass = isequal(optimum1,optimum2,optimum3);
 
maxerr = NaN;

if opt.Display
    plot(1:3,metric)
    set(gca,'xtick',[1 2 3],'xticklabel',{'Exponential','Stretched Exp.','Prod. Stretched Exp.'})
    xtickangle(45)
    ylabel('AICc')
    grid on, axis tight, box on
end

end

