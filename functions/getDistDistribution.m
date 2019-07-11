function Distribution = getDistDistribution(Signal,Method,Opts)

if nargin>3
  if ~isa(Opts,'DAoptions')
    error('Third argument must a valid DAoptions class object.')
  end
end

switch lower(Method)
    case 'apt'
        [Distribution,DistanceAxis] = APT(Signal.DipEvoFcn,Signal.TimeStep,Opts);
    case {'tikhonov','tv','huber'}
        [Kernel,DistanceAxis] = getKernel(Signal.Length,Signal.TimeStep,[],[],sqrt(Signal.Background));
        Distribution = regularize(Signal.ClusterFcn,Kernel,Method,Opts);
    case {'obir tikhonov','obir tv'}
        Method = Method(5:end);
        [Kernel,DistanceAxis] = getKernel(Signal.Length,Signal.TimeStep,[],[],sqrt(Signal.Background));
        Distribution = OBIR(Signal.ClusterFcn,Kernel,Method,Opts);
end

Distribution = DAdistribution('DistanceAxis',DistanceAxis,...
                          'Distribution',Distribution,...
                          'Signal',Signal.DipEvoFcn,...
                          'TimeAxis',Signal.TimeAxis,...
                          'SignalID',Signal.ID);