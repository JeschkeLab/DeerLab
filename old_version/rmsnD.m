function rms = rmsnD(v,x,y,logB,handles)
%RMS1D	Root mean square error of function exp(-k*t^(n/3)).
%	rms = rmsnD(v,x,y).
%	
%  Parameter: v(1) Amplitude, v(2) Zeitkonstante, v(3) Basislinie

%	Copyright (c) 1998 by Gunnar Jeschke

if v(1)<0, rms=1.0e6; return; end;
sim=decaynD(v,x,handles.hom_dim);
diff=sim-y;
rms=sum(diff.*diff);

		

