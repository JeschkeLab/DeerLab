function ampnoise = whitenoise(N,level,seed)

if nargin<3
    seed  = 2;
end
rng(seed);
noise = randn(N,1);
ampnoise = 0;
amp = 0;
while std(ampnoise)<level
    amp = amp+0.0001;
    ampnoise = amp*noise;
end

end