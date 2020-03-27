function [pass,maxerr] = test(opt)

% Check indifference of correctphase() towards input dimensionality

S = rand(1,100) + 1i*rand(1,100);
[Vr1,Vi1] = correctphase(S);
[Vr2,Vi2] = correctphase(S.');

% Pass 1: all real part signals are equal
pass(1) = isequal(Vr1,Vr2);
% Pass 2: all imaginary part signals are equal
pass(2) = isequal(Vi1,Vi2);
% Pass 3: all real part signals are columns
pass(3) = iscolumn(Vr1) | iscolumn(Vr2);
% Pass 3: all imaginary part signals are columns
pass(4) = iscolumn(Vi1) | iscolumn(Vi2);

pass = all(pass);

maxerr = NaN;
 

end