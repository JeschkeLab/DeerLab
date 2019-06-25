function B0 = get_resonance_field(freq,g)
% B0 = get_resonance_field(freq,g)
%
% freq  frequency [GHz]
% g     g value, optional, defaults to 2.0059 for nitroxide
% 
% B0    resonance field [G]

if ~exist('g','var')
    g = 2.0059; % nitroxide 
end

B0 = 1e9*1e4*hbar*2*pi*freq/(bmagn*g);
