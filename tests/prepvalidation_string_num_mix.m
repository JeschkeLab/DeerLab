function [err,data,maxerr] = test(data,opts)

Parameters = validationParameters();

Parameters(1).name = 'PenaltyType';
Parameters(1).vector = {'tikhonov','tv','huber'};

Parameters(2).name = 'RegParam';
Parameters(2).range = [10 100];
Parameters(2).trials = 10;

Parameters(3).name = 'NoiseLevel';
Parameters(3).range = [0.05 0.1];
Parameters(3).trials = 3;
Parameters(3).sampling = 'rand';

Parameters(4).name = 'SelectionMethod';
Parameters(4).vector = {'aic','gcv','gml'};

output = prepvalidation(Parameters);

err = ~isequal(size(output),[3*10*3*3 4]);
data = [];
maxerr = [];
end