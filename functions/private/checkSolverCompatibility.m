function checkSolverCompatibility(Solver,Method)
switch Method
    case 'tikhonov'
        allowedInput = {'fnnls','fmincon','cvx','lsqnonneg','bppnnls'};
    case {'tv','huber'}
        allowedInput = {'cvx','fmincon'};
    case 'custom'
        allowedInput = {'fmincon'};
end
if ~any(strcmp(allowedInput,Solver))
    error('The ''%s'' solver is not compatible with the ''%s'' method.',Solver,Method)
end
end