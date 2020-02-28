function [pass,maxerr] = test(opt)

% Check that fitbackground() can fit an exponential in log-scale

t = linspace(0,5,100).';
d = 3;

Bmodel = @(k)exp(-(k*t).^(d/3));

k = [0.5 1 1.5];

for ik = 1:numel(k)
  Bs{ik} = Bmodel(k(ik));
  B = Bs{ik};
  tstart = t(20);
  Bfit{ik} = fitbackground(B,t,@td_exp,tstart,'Logfit',true);
  maxerr(ik) = max(abs(Bfit{ik} - Bs{ik}));
end

maxerr = max(maxerr);

pass = maxerr < 1e-5;

%Plot results
if opt.Display
  figure,clf
  for ik = 1:numel(k)
    subplot(1,k,ik)
    plot(t,Bs{ik},t,Bfit{ik})
  end
end

end
