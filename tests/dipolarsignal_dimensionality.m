function [err,data,maxerr] = test(opt,olddata)

t = linspace(-1,4,20);
r = linspace(1,6,50);
B = td_exp(t,0.3);
lam = 0.5;
P = rd_onegaussian(r,[3 0.3]);

V1 = dipolarsignal(t,r,P,'moddepth',lam,'background',B);
V2 = dipolarsignal(t.',r,P,'moddepth',lam,'background',B);
V3 = dipolarsignal(t,r.',P,'moddepth',lam,'background',B);
V4 = dipolarsignal(t,r,P.','moddepth',lam,'background',B);
V5 = dipolarsignal(t.',r,P.','moddepth',lam,'background',B.');
V6 = dipolarsignal(t.',r.',P,'moddepth',lam,'background',B);
V7 = dipolarsignal(t.',r,P,'moddepth',lam,'background',B.');
V8 = dipolarsignal(t,r.',P.','moddepth',lam,'background',B);
V9 = dipolarsignal(t.',r.',P.','moddepth',lam,'background',B.');


err(1) = ~isequal(V1,V2,V3,V4,V5,V5,V6,V7,V8,V9);
err(2) = ~iscolumn(V1) | ~iscolumn(V2) | ~iscolumn(V3) | ~iscolumn(V4) | ~iscolumn(V5) | ~iscolumn(V6) | ~iscolumn(V7) | ~iscolumn(V8) | ~iscolumn(V9);

err = any(err);
maxerr = max(max(abs(V1 - V2)));
data = [];

end