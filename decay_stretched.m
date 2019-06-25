function sim=decay_stretched(v,x);
%
%
sim=v(2)*exp(-v(1)*x.^(v(3)/3));
