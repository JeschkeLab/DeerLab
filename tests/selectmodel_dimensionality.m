function [err,data,maxerr] = test(opt,olddata)

t = linspace(0,5,80);
S = dipolarsignal(t,3);
r = linspace(1,6,50);
K = dipolarkernel(t,r);

[opt1,f1,params1] = selectmodel({@rd_onegaussian,@rd_fourgaussian},S,r,K,'aicc');
[opt2,f2,params2] = selectmodel({@rd_onegaussian,@rd_fourgaussian},S.',r,K,'aicc');
[opt3,f3,params3] = selectmodel({@rd_onegaussian,@rd_fourgaussian},S,r.',K,'aicc');
[opt4,f4,params4] = selectmodel({@rd_onegaussian,@rd_fourgaussian},S.',r.',K,'aicc');

err(1) = ~iscolumn(f1) | ~iscolumn(f2) | ~iscolumn(f3) | ~iscolumn(f4);
err(2) = iscolumn(params1) | iscolumn(params2) | iscolumn(params3)| iscolumn(params4);
err(2) = ~isequal(opt1,opt2,opt3,opt4);

err = any(err);
maxerr = max(opt1 - opt2);
data = [];

end