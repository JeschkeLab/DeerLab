function [err,data,maxerr] = test(opt,olddata)


r = linspace(0,50,500);

%Control row/column effect
P1 = rd_tworice(r  ,[2.5 0.4 4.0 0.4 0.5]);
P2 = rd_tworice(r.',[2.5 0.4 4.0 0.4 0.5]);
err(1) = ~isequal(P1,P2);

%Control non-negativity
err(2) = any(P1<0);

P3 = rd_tworice(r,[1.0 0.1 1.0 0.1 0.5]);
P4 = rd_tworice(r,[10 5 10 5 0.5]);

%Control non-negativity ov default boundaries
err(3) = any(P3<0) | any(P4<0);

%Control there are no NaN
err(4) = any (isnan(P1)) | any(isnan(P2)) | any(isnan(P3)) | any(isnan(P4));

err = any(err);
data = [];
maxerr = NaN;

end