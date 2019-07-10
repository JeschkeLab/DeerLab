function [err,data] = test(opt,olddata)

%======================================================
% Test daopts class rects well to inputs
%======================================================

err = false;

try
DAoptions('BackgroundModel','exponential');
err(1) = false;
catch
err(1) = true;    
end

try
DAoptions('BackgroundModel','wronginput');
err(2) = true;
catch
err(2) = false;
end


err = any(err);
data = [];

end