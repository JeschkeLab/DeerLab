function [pass,maxerr] = test(opt)

% Check that sensitivan behaves correctly with arrays of changing size

rng(1)
param.N = [300 250 500];

% Pass: an error is returned as expected
try
    sensitivan(@(param)myfcn(param),param);
    pass = false;
catch
    pass = true;
end

maxerr = NaN;

end


function [x] = myfcn(p)

N = p.N;

rng(1)
x = rand(N,1);

end