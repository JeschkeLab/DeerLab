function [err,data,maxerr] = test(opt,olddata)

varpar.x = linspace(1,8,3);
[stats1] = sensitivan(@fcn,varpar);
varpar.x = varpar.x.';
[stats2] = sensitivan(@fcn,varpar);

err(1) = ~iscolumn(stats1(1).median) | ~iscolumn(stats2(1).median) | ~iscolumn(stats1(2).median) | ~iscolumn(stats2(2).median);
err(2) = ~isequal(stats1,stats2);

err = any(err);
maxerr = max(stats1(1).median - stats2(1).median);
data = [];


    function [y1,y2] = fcn(varpar)
        rng(varpar.x)
        y1 = linspace(0,5,20);
        rng(varpar.x)
        y2 = linspace(0,5,20);
    end

end

