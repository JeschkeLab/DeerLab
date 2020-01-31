function [err,data,maxerr] = test(opt,olddata)


r = linspace(0,50,500);

%Control row/column effect
P1 = rd_onegaussian(r,[3,0.5]);
P2 = rd_onegaussian(r.',[3,0.5]);
err(1) = ~isequal(P1,P2);

%Control non-negativity
err(2) = any(P1<0);

P3 = rd_onegaussian(r.',[1,0.2]);
P4 = rd_onegaussian(r.',[20,5]);

%Control non-negativity ov default boundaries
err(3) = any(P1<0) | any(P2<0);

%Control there are no NaN
err(4) = any (isnan(P1)) | any(isnan(P2)) | any(isnan(P3)) | any(isnan(P4));

err = any(err);
data = [];
maxerr = NaN;

end