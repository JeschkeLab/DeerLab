function [pass,maxerr] = test(opt)

param.N = [5 5 7];
param.M = [7 6 2];

rng(1)
try
    [mean]  = sensitivan(@(param)myfcn(param),param);
    err = true;
catch
    err = false;
end

 
maxerr = 0;

end


function [x] = myfcn(p)

N = p.N;
M = p.M;

rng(1)
x = rand(N,M);



end