function [pass,maxerr] = test(opt)

% Check that regoperator returns correct sizes

r = 1:30;
n = numel(r);
L0 = regoperator(r,0);
L1 = regoperator(r,1);
L2 = regoperator(r,2);

% Pass: all operators have the right sizes
sizeok = @(M,sz)all(size(M)==sz);
pass = sizeok(L0,[n n]) && sizeok(L1,[n-1 n]) && sizeok(L2,[n-2 n]);

maxerr = NaN;

end