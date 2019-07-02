function rms = rms_user_fct(v,x,y,userfunction)
%RMS2D	Root mean square error of user background function fct(k,t)
%	rms = rms_user_fct(v,x,y,fct).
%	
%  Parameters: v(2) Amplitude, v(1) decay constant k
%  fct is a Matlab expression that must contain variable k and t and must
%  not contain any other variables, amplitude scaling is external to this
%  function

%	Copyright (c) 2004 by Gunnar Jeschke

if v(1)<0, rms=1.0e6; return; end;
sim=decay_user_fct(v,x,userfunction);
diff=sim-y;
rms=sum(diff.*diff);

		

