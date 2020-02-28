function [pass,maxerr] = test(opt)


r = linspace(0,50,500);

%Control row/column effect
P1 = rd_wormgauss(r,[3.7 10 0.2]);
P2 = rd_wormgauss(r.',[3.7 10 0.2]);
err(1) = ~isequal(P1,P2);

%Control non-negativity
err(2) = any(P1<0);

P3 = rd_wormgauss(r,[1.5 2 0.001]);
P4 = rd_wormgauss(r,[10 100 2]);

%Control non-negativity ov default boundaries
err(3) = any(P1<0) | any(P2<0);

%Control there are no NaN
err(4) = any (isnan(P1)) | any(isnan(P2)) | any(isnan(P3)) | any(isnan(P4));


pass = all(err);
 
maxerr = NaN;

end