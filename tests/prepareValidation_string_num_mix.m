function [err,data,maxerr] = test(data,opts)

Parameters = validationParameters();

Parameters(1).name = 'PenaltyType';
Parameters(1).vector = {'tikhonov','tv','huber'};

Parameters(2).name = 'RegParam';
Parameters(2).range = [10 100];
Parameters(2).trials = 10;

Parameters(3).name = 'RegMatrixOrder';
Parameters(3).vector = [1 2 3];
Parameters(3).trials = 3;

Parameters(4).name = 'Solver';
Parameters(4).vector = {'fnnls','fmincon','lsqnonneg'};
Parameters(4).trials = 3;

output = prepareValidation(Parameters);

err = ~isequal(size(output),[3*10*3*3 4]);
data = [];
maxerr = [];
end