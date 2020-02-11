function [err,data,maxerr] = test(opt,olddata)


r = linspace(0,50,500);

%Control row/column effect
P1 = rd_spheresurf(r,3.5);
P2 = rd_spheresurf(r.',3.5);
err(1) = ~isequal(P1,P2);

%Control non-negativity
err(2) = any(P1<0);

P3 = rd_spheresurf(r.',0.5);
P4 = rd_spheresurf(r.',30);

%Control non-negativity of default boundaries
err(3) = any(P1<0) | any(P2<0);

%Control there are no NaN
err(4) = any (isnan(P1)) | any(isnan(P2)) | any(isnan(P3)) | any(isnan(P4));

%Control non-negativity
dr = mean(diff(r));

%Control normalization
err(5) = round(sum(P1)*dr)~=1 |round(sum(P2)*dr)~=1 | round(sum(P3)*dr)~=1 | round(sum(P4)*dr)~=1;

err = any(err);
data = [];
maxerr = NaN;

end