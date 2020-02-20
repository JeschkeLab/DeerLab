function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Test RIDME overtone coefficients
%=======================================
n = 5;
Tmix = 15; % us
T1 = 88; % us

TheoreticalCoeff = 1-exp(-Tmix/T1./(1:n));
TheoreticalCoeff = TheoreticalCoeff/sum(TheoreticalCoeff);
CalculatedCoeff = overtones(n,Tmix,T1);
error1 = abs(CalculatedCoeff-TheoreticalCoeff);

err = any(error1 > 1e-10);

err = any(err);
maxerr = max(error1);
data = [];

end
