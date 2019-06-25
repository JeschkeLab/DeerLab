function sim=decay_user_fct(v,x,userfunction);
%
%
t=x;
k=v(1);
eval(userfunction);
sim=v(2)*sim;
