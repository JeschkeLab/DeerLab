function [pass,maxerr] = test(opt)

% Check indifference of sensitivan() towards input dimensionality

varpar.x = linspace(1,8,3);
[stats1] = sensitivan(@fcn,varpar);
varpar.x = varpar.x.';
[stats2] = sensitivan(@fcn,varpar);

% Pass 1: all medians are columns
pass(1) = iscolumn(stats1(1).median) & iscolumn(stats2(1).median) & iscolumn(stats1(2).median) & iscolumn(stats2(2).median);
% Pass 1: all statistics are equal
pass(2) = isequal(stats1,stats2);

pass = all(pass);

maxerr = NaN;


    function [y1,y2] = fcn(varpar)
        rng(varpar.x)
        y1 = linspace(0,5,20);
        rng(varpar.x)
        y2 = linspace(0,5,20);
    end

end

