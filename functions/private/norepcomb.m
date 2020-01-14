
function output = norepcomb(vals,taken)

[~,Nvals] = size(vals);

Nrows = 2.^(Nvals);
Ncycles = Nrows;

for i = 1:Nvals
    settings = (0:1);
    Ncycles = Ncycles/2;
    nreps = Nrows./(2*Ncycles);
    settings = settings(ones(1,nreps),:);
    settings = settings(:);
    settings = settings(:,ones(1,Ncycles));
    Indexes(:,Nvals-i+1) = settings(:);
end

Index = Indexes(sum(Indexes,2) == taken,:);
nrows = size(Index,1);
[Nrows,~] = find(Index');
output = reshape(vals(Nrows),taken,nrows).';

end
