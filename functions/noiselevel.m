function level = noiselevel(Signal)

N = length(Signal);
cutoff = Signal(ceil(4/5*N):N);
cutoff = cutoff - mean(cutoff);
level = std(cutoff);

end