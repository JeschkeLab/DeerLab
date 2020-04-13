%
% GETMAXARGOUT Determine number of outputs returned by function handle
%
%   varargout = GETMAXARGOUT(fcnHandle,argin)
%   Determines the number of outputs which the function handle (fcnHandle)
%   is gonna return when passing the inputs in (argin).
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function nout = getmaxnargout(fcnHandle,argin)
nout = nargout(fcnHandle);
variableOutputs = nout<0;

if ~variableOutputs
    return;
end

% If the functions defines a variable number of outputs, iteratively increase
% number of requested outputs until the function crashes, to determine maximum
% number of outputs.
nout = abs(nout)-1;
done = false;
while ~done
    try
        nout = nout+1;
        varargout = cell(1,nout);
        [varargout{:}] = fcnHandle(argin);
    catch
        nout = nout-1;
        done = true;
    end
end
end