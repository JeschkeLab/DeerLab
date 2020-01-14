%
% TD_DMPDEER Time-domain dipolar multi-pathway (DMP) model for multi-pulse DEER
%
%       V = TD_DMPDEER(t,r,P,taus,ts,p)
%       Computes the dipolar signal (given a distance distribution (P) on a
%       distance axis (r)) originating from any multi-pulse DEER
%       sequence with (etas) and (ts) pulse sequence parameters
%
%               Label = 1               Label = 2               Label = N
%        <----------------------><------------------------>   ...-->    
%         _         _      _               _      _
%        | |  tau1 | | t1 |||tau1-t1      | |    |||                 /\
%        | |<----->| |<-->|||<--><------->| |<-->|||<-----><--...-->/  \
%        | |       | |    |||       tau2  | | t2 ||| tau2-t2       /    \
%     ---------------------------------------------------------------------
%       where taus = [tau1 tau2 ... tauN] and ts = [t1 t2 ... tN]. Each pump 
%       pulse has a probability (p) to invert the magnetization.       
%
%       V = TD_DMPDEER(t,r,P,taus,ts,p,Bpar)
%       In order to simulate the signal with background, the decay rate (k)
%       and dimensionality (d) of a stretched expnential function can be
%       passed as an array (Bpar), where Bpar=[k d].
%
%       V = TD_DMPDEER(t,r,P,taus,ts,p,Bpar,lab)
%       While each pump-pulse gets a label (1,2,...,N) assigned according to the
%       tau-block it resides in, the labels of the pump pulses can be
%       specified by means of a last argument (lab). This allows, e.g. to
%       simulate where the observer pulses have a finite probability to
%       invert the magnetization and hence to act as pump pulses.
%

% This file is a part of DeerAnalysis. License is MIT (see LICENSE.md).
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function V = td_dmpdeer(t,r,P,taus,ts,prob,Bparam,labels)

if nargin<6
   error('Not enough input arguments.') 
end
if ~iscell(ts)
    validateattributes(ts,{'numeric'},{'nonempty'},mfilename,'taus')
    ts = mat2cell(ts,numel(ts),1);
else
    for i=1:length(ts)
        if ~iscolumn(ts{i})
            ts{i} = ts{i}.';
        end
    end
end
%Input validation
validateattributes(t,{'numeric'},{'nonempty'},mfilename,'t')
validateattributes(r,{'numeric'},{'nonnegative','nonempty'},mfilename,'r')
validateattributes(P,{'numeric'},{'nonnegative','nonempty'},mfilename,'P')
validateattributes(taus,{'numeric'},{'nonempty'},mfilename,'taus')
validateattributes(prob,{'numeric'},{'nonempty'},mfilename,'prob')

if nargin<7 || isempty(Bparam)
    Bparam = [0 3];
else
    validateattributes(Bparam,{'numeric'},{'nonnegative'},mfilename,'Bparam')
end

if nargin<8
    labels = [];
else
    validateattributes(labels,{'numeric'},{'nonnegative','nonempty'},mfilename,'labels')
end

if numel(taus)~=numel(ts)
    error('The number of tau-times and t-times must be equal.')
end
if numel(prob)>1 && numel(prob)~=numel(ts)
    error('The number inversion probabilities must be equal to the number of taus.')
end
if ~iscell(ts)
    ts = mat2cell(ts);
end
if numel(prob)==1
   prob = prob*ones(numel(taus),1);
end
if ~iscolumn(t)
    t = t.';
end

%Get number of pump pulses used
Npumps = numel(taus);

%Get pump pulse labels if not given
if isempty(labels)
    labels = 1:Npumps;
end

%Prepare containers
V = zeros(numel(t),1);
B = ones(numel(t),1);

%Get background parameters
k = Bparam(1);
d = Bparam(2);

%Probability of not inverting with any pulse
lam0  = prod(1-prob);
V = V + lam0;

%Loop over number of pump pulses inverting the magnetization
for Npumped = 1:Npumps
    
    %Get all possible (non-repeated) combinations
    combs = norepcomb(1:Npumps,Npumped);
    
    %Loop over pathways involving Npumped-inversions
    for n=1:size(combs,1)
        
        %Get current pathway label
        comb = combs(n,:);
        notcomb = setdiff(1:Npumps,comb);

        %Initialize variables
        s = 1;
        prev = 0;
        tevo = 0;
       
        %Probability of current dipolar pathway
        if ~isempty(notcomb)
            lambda =  prod(prob(comb))*prod(1-prob(notcomb));
        else
            lambda =  prod(prob(comb));
        end
        
        %Loop over the events in the pathway
        for j=1:numel(comb)
            %Compute current phase sign
            s = s*(-1)^(1 + labels(comb(j)) - prev);
            prev = labels(comb(j));
            %Accumulate the dipolar evolution time for the pathway
            tevo = tevo + s*(taus(comb(j)) - ts{comb(j)});
        end
        
        %Just some formatting of the vectors for the dipolarsignal function
        if mean(diff(tevo))<0
            tevo =  -tevo;
        end
        
        %Add the dipolar pathway contribution to the signal
        V = V + lambda*dipolarsignal(tevo,r,P);
        
        %Add the dipolar pathway contribution to the background
        B = B.*exp(-(k*lambda*abs(tevo)).^(d/3));
    end
    
end

    %Apply the background
    V = V.*B;

    
    return
    
    
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

return