function [err,data,maxerr] = test(opt,oldata)

Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam  = [3,0.5];
P = rd_onegaussian(r,InputParam);

K = dipolarkernel(t,r);
S = K*P;

InitialGuess = [2 0.1];
%nm
Pfit1 = fitparamodel(S,@rd_onegaussian,r,K,InitialGuess,'Solver','fmincon');
%A
r = r*10;
Pfit2 = fitparamodel(S,@rd_onegaussian,r,K,InitialGuess,'Solver','fmincon');

err = any(abs(Pfit1 - Pfit2)>1e-12);
maxerr = max(abs(Pfit1 - Pfit2));
data = [];

if opt.Display
   figure(1),clf,hold on
   plot(t,Pfit1,'b')
   plot(t,Pfit2,'r')
end

end