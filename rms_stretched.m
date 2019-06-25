function rms = rms_stretched(v,x,y,handles)
%RMS_streched	Root mean square error of function exp(-k*t^ksi).
%	rms = rms_stretched(v,x,y).
%	
%  Parameter: v(2) Amplitude, v(1) Zeitkonstante, v(3) ksi

%	Copyright (c) 2004 by Gunnar Jeschke

if v(1)<0, rms=1.0e6; return; end;

% distr=handles.Pake_r.^(v(3)-1);
% distr=0.01*distr/sum(distr);
% logB=distr*handles.Pake_kernel;
sim=decaynD(v(1:2),x,v(3));
% sim=decay_stretched(v,x);
diff=sim-y;
rms=sum(diff.*diff);

		

