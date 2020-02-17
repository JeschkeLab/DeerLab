function [err,data,maxerr] = test(opt,olddata)

t = linspace(0,5,30);
d = 3;
k = 0.5;
bckg = exp(-(k*t).^(d/3));

clear aptkernel
t1 = backgroundstart(bckg,t,@td_exp);
clear aptkernel
t2 = backgroundstart(bckg.',t,@td_exp);
clear aptkernel
t3 = backgroundstart(bckg.',t.',@td_exp);


err(1) = ~isequal(t1,t2,t3);
err(2) = ~iscolumn(t1) | ~iscolumn(t2) | ~iscolumn(t3);

err = any(err);

maxerr = max(abs(t1 - t2));
data = [];

end