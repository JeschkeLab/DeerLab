function [pass,maxerr] = test(opt)

% Compare td_dmpdeer output to explicit calculation using dipolarkernel

t = linspace(-0.3,3,401);
r = linspace(1,3,201);

r0 = 2;
w = 0.2;
P = rd_onegaussian(r,[r0 w]);

pathinfo = [0.5 0; 0.2 1.2];

Bfun = @td_exp;

V1 = td_dmpdeer(t,r,P,pathinfo,Bfun);
K = dipolarkernel(t,r,pathinfo,Bfun);
V2 = K*P(:);

maxerr = max(abs(V1-V2));
pass = maxerr < 1e-10;
 

if opt.Display
  plot(t,V1,t,V2);
  legend('td_dmpdeer','dipolarkernel');
end

end
