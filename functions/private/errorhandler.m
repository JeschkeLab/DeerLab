%
% ERRRORHANDLER Wrapper for formating DeerLab exceptions
%
%   out = ERRRORHANDLER(fcnh,env,inputs)
% 

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function output = errorhandler(fcnh,env,varargin)
   try
     % Evaluate the function handle  
     output = feval(fcnh, varargin{:});
     
   % If there is an errorcatch it and format it properly
   catch err
       
     % Check where the error has ocurred
     switch env
         case 'modelfcn'
             msg = sprintf('An error was found in the model function. \n\nCaused by:\n  %s',err.message);
     end
     
     % Construct error
     newerr.message = msg;
     newerr.stack(1,1) = err.stack(1);
     newerr.stack(2,1) = err.stack(end-1);
     error(newerr)
   end
end