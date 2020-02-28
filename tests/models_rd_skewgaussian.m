function [pass,maxerr] = test(opt)


r = linspace(0,50,500);

%Control row/column effect
P1 = rd_skewgaussian(r,[3,0.5,5]);
P2 = rd_skewgaussian(r.',[3,0.5,5]);
err(1) = ~isequal(P1,P2);

%Control non-negativity
err(2) = any(P1<0);

P3 = rd_skewgaussian(r.',[1,0.2,-15]);
P4 = rd_skewgaussian(r.',[20,5,15]);

%Control non-negativity ov default boundaries
err(3) = any(P1<0) | any(P2<0);

%Control there are no NaN
err(4) = any (isnan(P1)) | any(isnan(P2)) | any(isnan(P3)) | any(isnan(P4));

pass = all(err);
 
maxerr = NaN;

end