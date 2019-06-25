function rms = rms_ubckg(v,y,bckgdat)
%RMS1D	Root mean square error of scaled background data from data y.
%	rms = rms1D(v,x,y).
%	
%  Parameter: v(1) Amplitude, v(2) Zeitkonstante, v(3) Basislinie

%	Copyright (c) 1998 by Gunnar Jeschke

if v(1)<0, rms=1.0e6; return; end;
sim=v(2)*exp(v(1)*log(bckgdat));
diff=sim-y;
rms=sum(diff.*diff);

		

