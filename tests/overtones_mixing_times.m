function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Test RIDME overtone coefficients
%=======================================
Order = 5;

Tmix = 15; %us
T1 = 88; %us
TheoreticalCoeff = [0.42,0.22,0.15,0.11,0.09];
CalculatedCoeff = overtones(Order,Tmix,T1);
error1 = abs(TheoreticalCoeff - CalculatedCoeff);


Tmix = 50; %us
T1 = 88; %us
TheoreticalCoeff = [0.40,0.23,0.16,0.12,0.10];
CalculatedCoeff = overtones(Order,Tmix,T1);
error2 = abs(TheoreticalCoeff - CalculatedCoeff);


Tmix = 500; %us
T1 = 88; %us
TheoreticalCoeff = [0.24,0.22,0.20,0.18,0.16];
CalculatedCoeff = overtones(Order,Tmix,T1);
error3 = abs(TheoreticalCoeff - CalculatedCoeff);


err(1) = any(error1 > 1e-2);
err(2) = any(error2 > 1e-2);
err(3) = any(error3 > 1e-2);
err = any(err);
maxerr = max(error1);
data = [];

end