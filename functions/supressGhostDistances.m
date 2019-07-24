function Signal = supressGhostDistances(Signal,NRadicals)

if nargin<3
    NRadicals = 2;
end

Scaling = 1/(NRadicals - 1);

Signal = Signal.^Scaling;

end





