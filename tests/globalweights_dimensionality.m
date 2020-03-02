function [pass,maxerr] = test(opt)

% Check indifference of globalweights() towards input dimensionality

S1 = ones(40,1);
S2 = ones(70,1);
S3 = ones(120,1);
S = {S1,S2,S3};
levels = [0.1 0.05 0.09];

w1 = globalweights(S,levels);
w2 = globalweights(S.',levels);
w3 = globalweights(S,levels.');
w4 = globalweights(S.',levels.');

% Pass 1: all weights are equal
pass(1) = isequal(w1,w2,w3,w4); 
% Pass 2: all weights are column vectors
pass(2) = isequal(w1,w2,w3,w4); 

pass = all(pass);

maxerr = NaN;

end
