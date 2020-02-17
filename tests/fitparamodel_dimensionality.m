function [err,data,maxerr] = test(opt,oldata)

t = linspace(0,5,50);
r = linspace(1,6,50);
P = rd_onegaussian(r,[4 0.3]);

K = dipolarkernel(t,r);
S = K*P;

[~,Pfit1] = fitparamodel(S,@rd_onegaussian,r,K);
[~,Pfit2] = fitparamodel(S.',@rd_onegaussian,r,K);
[~,Pfit3] = fitparamodel(S,@rd_onegaussian,r.',K);
[~,Pfit4] = fitparamodel(S.',@rd_onegaussian,r.',K);


mymodel = @(t,param)K*rd_onegaussian(r,param);

param0 = [4 0.3];

[~,Vfit1] = fitparamodel(S,mymodel,t,param0,'upper',[100 100],'lower',[0 0]);
[~,Vfit2] = fitparamodel(S.',mymodel,t,param0.','upper',[100 100],'lower',[0 0]);
[~,Vfit3] = fitparamodel(S,mymodel,t.',param0,'upper',[100 100],'lower',[0 0]);
[~,Vfit4] = fitparamodel(S.',mymodel,t.',param0.','upper',[100 100],'lower',[0 0]);

err(1) = ~isequal(Pfit1,Pfit2,Pfit3,Pfit4);
err(2) = ~iscolumn(Pfit1) | ~iscolumn(Pfit2) | ~iscolumn(Pfit3) | ~iscolumn(Pfit4);
err(3) = ~iscolumn(Vfit1) | ~iscolumn(Vfit2) | ~iscolumn(Vfit3) | ~iscolumn(Vfit4);

err = any(err);

maxerr = max(abs(Pfit1 - Pfit2));
data = [];


end