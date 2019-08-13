function weights = globalweights(Signal)

N = cellfun('length',Signal);
NoiseLevel = cellfun(@noiselevel,Signal);
weights = sum(N.*NoiseLevel)./(N.*NoiseLevel);
weights = weights/sum(weights);

end