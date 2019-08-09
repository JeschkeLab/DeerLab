function ampnoise = whitenoise(N,level,seed)

if nargin<3
    seed  = 2;
end
rng(seed);
noise = rand(N,1);
noise = noise - mean(noise);
noise = noise/max(noise);
amp = 0;
ampnoise = 0;
while std(ampnoise)<level
    amp = amp+0.001;
    ampnoise = amp*noise;
end

end