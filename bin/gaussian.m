function y = gaussian(x,x0,w)
    y = exp(-(x  - x0).^2/w.^2);
end