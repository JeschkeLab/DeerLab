function rms=rms_phi(phi,tr);
% r.m.s. of imaginary part after phase correction
itr=imag(tr*exp(i*phi));
rms=sqrt(sum(itr.*itr));
