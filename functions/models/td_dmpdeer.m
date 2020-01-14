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
    Bparam = [0 1];
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
        B = B.*exp(-(k*lambda*abs(tevo)).^d);
        
    end
    
    %Apply the background
    V = V.*B;
    
end
