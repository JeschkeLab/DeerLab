function [err,data,maxerr] = test(opt,olddata)

t = linspace(0,5,100);
d = 3;
k = 0.5;
bckg = exp(-(k*t).^(d/3));

tstart = [t(1) t(end)];

[B1,~,param1,~] = fitbackground(bckg,t,@td_strexp,tstart);
[B2,~,param2,~] = fitbackground(bckg.',t,@td_strexp,tstart);
[B3,~,param3,~] = fitbackground(bckg,t.',@td_strexp,tstart);
[B4,~,param4,~] = fitbackground(bckg.',t.',@td_strexp,tstart.');


err(1) = ~isequal(B1,B2,B3,B4);
err(2) = ~iscolumn(B1) | ~iscolumn(B2) | ~iscolumn(B3) | ~iscolumn(B4);
err(3) = iscolumn(param1) | iscolumn(param2) | iscolumn(param3) | iscolumn(param4);

err = any(err);

maxerr = max(abs(B1 - B2));
data = [];

end