function [pass,maxerr] = test(opt)

param.N = [300 250 500];

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

rng(1)
x = rand(N,1);

end