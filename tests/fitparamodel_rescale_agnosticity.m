function [pass,maxerr] = test(opt)

% Check that fitparamodel rescaling does not change the results

rng(1)

t = linspace(0,5,100);
r = linspace(2,6,100);
P = dd_gauss(r,[4 0.5]);
B = bg_exp(t,0.3);
K = dipolarkernel(t,r,0.3,B);

scale = 1e9;
V = K*P + whitegaussnoise(t,0.005);
InitialGuess = [0.5 0.5 2 0.3];

lower = [0 0 1 0.1];
upper = [1 1 20 5];
mymodel = @(t,param)dipolarkernel(t,r,param(1),bg_exp(t,param(2)))*dd_gauss(r,param(3:4));

parfit1 = fitparamodel(V*scale,mymodel,t,InitialGuess,'Rescale',true,'MultiStart',5,'Lower',lower,'Upper',upper);
parfit2 = fitparamodel(V,mymodel,t,InitialGuess,'Rescale',false,'MultiStart',5,'Lower',lower,'Upper',upper);

Pfit1 = dd_gauss(r,parfit1(3:4));
Pfit2 = dd_gauss(r,parfit2(3:4));

%Pass 1: distance distribution is well fitted
pass(1) = all(abs(Pfit1 - P) < 2e-1);
pass(2) = all(abs(Pfit2 - Pfit1) < 1e-2);

pass = all(pass);

maxerr = max(abs(Pfit1 - P));

if opt.Display
   plot(r,P,'k',r,Pfit1,'r',r,Pfit2,'b')
   legend('truth','rescaled','normalized')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end