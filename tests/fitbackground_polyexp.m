function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Exponential background fit
%======================================================

t = linspace(0,5,100).';
d = 3;

bg = @(k)exp(-(k*t).^(d/3));

k = [0.5 1 1.5];

for ik = 1:numel(k)
  bckg{ik} = bg(k(ik));
  data2fit = bckg{ik};
  tstart = t(20);
  fit{ik} = fitbackground(data2fit,t,@td_exp,tstart,'Logfit',true);
  maxerr(ik) = max(abs(fit{ik} - bckg{ik}));
end
maxerr = max(maxerr);
err = maxerr>1e-5;

data = [];

if opt.Display
  figure,clf
  for ik = 1:numel(k)
    subplot(1,k,ik)
    plot(t,bckg{ik},t,fit{ik})
  end
end

end
