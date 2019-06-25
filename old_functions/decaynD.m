function sim = decaynD(v,x,hom_dim)

sim = v(2)*exp(-abs(v(1)*x).^(hom_dim/3));
