function [pass,maxerr] = test(opt)

n = 5;
Tmix = 15; % us
T1 = 88; % us

TheoreticalCoeff = 1-exp(-Tmix/T1./(1:n));
TheoreticalCoeff = TheoreticalCoeff/sum(TheoreticalCoeff);
CalculatedCoeff = overtones(n,Tmix,T1);
error1 = abs(CalculatedCoeff-TheoreticalCoeff);

pass = all(error1 < 1e-10);

maxerr = max(error1);
 

end
