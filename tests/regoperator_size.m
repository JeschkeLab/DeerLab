function [pass,maxerr] = test(opt)

% Check that regoperator returns correct sizes

n = 30;
L0 = regoperator(n,0);
L1 = regoperator(n,1);
L2 = regoperator(n,2);

% Pass: all operators have the right sizes
sizeok = @(M,sz)all(size(M)==sz);
pass = sizeok(L0,[n n]) && sizeok(L1,[n-1 n]) && sizeok(L2,[n-2 n]);

maxerr = NaN;

end