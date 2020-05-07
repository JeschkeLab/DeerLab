function [pass,maxerr] = test(opt)

% Check that fitsignal() works with a dipolar evolution function

rng(1)

t = linspace(0,5,100);
r = linspace(2,6,30);
P = dd_gauss(r,[4.5 0.6]);
V = dipolarsignal(t,r,P,'noiselevel',0.01);

[Vfit1,Pfit1] = fitsignal(V,t,r,'P','none','none');
[Vfit2,Pfit2] = fitsignal(V.',t,r,'P','none','none');
[Vfit3,Pfit3] = fitsignal(V,t.',r,'P','none','none');
[Vfit4,Pfit4] = fitsignal(V,t,r.','P','none','none');
[Vfit5,Pfit5] = fitsignal(V.',t,r.','P','none','none');
[Vfit6,Pfit6] = fitsignal(V.',t.',r,'P','none','none');
[Vfit7,Pfit7] = fitsignal(V.',t.',r.','P','none','none');

% Pass 1: all signals are equal
pass(1) = isequal(Vfit1,Vfit2,Vfit3,Vfit4,Vfit5,Vfit6,Vfit7);
% Pass 2: all signals are columns
pass(2) = iscolumn(Vfit1) & iscolumn(Vfit2) & iscolumn(Vfit3) & iscolumn(Vfit4) & iscolumn(Vfit5) & iscolumn(Vfit6) & iscolumn(Vfit7);
% Pass 3: all distribution are equal
pass(3) = isequal(Pfit1,Pfit2,Pfit3,Pfit4,Pfit5,Pfit6,Pfit7);
% Pass 4: all distribution are columns
pass(4) = iscolumn(Pfit1) & iscolumn(Pfit2) & iscolumn(Pfit3) & iscolumn(Pfit4) & iscolumn(Pfit5) & iscolumn(Pfit6) & iscolumn(Pfit7);

pass = all(pass);

maxerr = NaN;

end