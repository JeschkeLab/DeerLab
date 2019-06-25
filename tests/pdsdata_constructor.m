function [err,data] = test(opt,olddata)


[t,data] = eprload('./examples/CT_DEER_5nm');
data = real(data);

formfactor = data/data(1);

trace = pdsdata();

trace.TimeAxis = t;
trace.RawData = data;

err = trace.Length~=755;
err = trace.TimeStep~=8;

data = [];