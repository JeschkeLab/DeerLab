%
% VALIDATE Statistical validation of results
%
%   [mean,std] = VALIDATE(p,vp)
%   Validates the parameter (p) according to the validation parameter
%   defintion in the structure (vp). The validation results are returned in
%   the form of the mean value (mean) and standard deviation (std) of the
%   parameter (p). The validation is run based on the code from the
%   script/function from which VALIDATE() has been called.
%
%   [mean,std] = VALIDATE(p,vp,filename)
%   If the script/function containing the code and parameters to be validated
%   are contained in an external file, this file name can be specified as a
%   third argument (filename);
%
%   [mean,std] = VALIDATE(p,vp,'Property',Value)
%   [mean,std] = VALIDATE(p,vp,filename,'Property',Value)
%   Additional (optional) arguments can be passed as property-value pairs.
%
% The properties to be passed as options can be set in any order.
%
%   'AxisHandle' - Axis handle to plot the state of the validation results
%                  at each parameter permutation (default = empty)
%
%   'RandPerm' - Specifies whether to randomly permutate the validation
%   `            parameters combinations (default = true)
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.


function [meanOut,stdOut] = validate(OutParam,Parameters,Filename,varargin)

if nargin<2
   error('Not enough input arguments. At least two input arguments required.') 
end

%Validate the required input attributes
validateattributes(OutParam,{'char'},{'nonempty'},mfilename,'OutParam')
validateattributes(Parameters,{'struct'},{'nonempty'},mfilename,'Parameters')

%Parse & validate required input
if ~isfield(Parameters,{'name','values'})
   error('The input structure must contain the ''name'' and ''values'' fields.') 
end

%-----------------------------------------------------
% Input parsening and validation
%-----------------------------------------------------

%Parse the case filename is given vs. not given
fileok = false;
if nargin>2
    [~,file] = fileparts(Filename);
    file = [file '.m'];
    if ~exist(file,'file')
        varargin = [{Filename} varargin];
    else
       fileok = true; 
    end
    scriptName = file;
end

%Parse optional inputs
[AxisHandle,RandPerm] = parseoptional({'AxisHandle','RandPerm'},varargin);

if isempty(RandPerm)
    RandPerm = true;
else
    validateattributes(RandPerm,{'logical'},{'nonempty'},mfilename,'RandPerm');
end

%If the script/function file has not been specified, then get the caller
if ~fileok
    %Get the call stack
    callstack = dbstack;
    filescalled = {callstack.file};
    %Get the name of the script/functions calling validate()
    scriptName = filescalled{2};
end

%If called via Run Section the passed live script cannot be opened
if contains(scriptName,'LiveEditorEvaluation')
    error('Validation cannot be called by ''Run Section''.');
end

%-----------------------------------------------------
% Preparation
%-----------------------------------------------------

%Prepare validation parameter permutations
validationParam = prepvalidation(Parameters,'RandPerm',RandPerm);

%-----------------------------------------------------
% Code scanning
%-----------------------------------------------------

%Open script or function
ID = fopen(scriptName);
%Prepare containers
code = {};
callerisfunction = false;
%Read through the code
while true
    useCode = true;
    %Read the next line line
    String = fgets(ID);
    %If end of the file reached, then stop reading
    if all(String==(-1)) || contains(String,'validate(')
        break;
    end
    %Check if the caller is a function...
    if contains(String,'function')
        %... if so get the names of the arguments being passed there
        callerisfunction = true;
        Pstart = findstr(String,'(');
        Pend = findstr(String,')');
        fcnInName = strsplit(String(Pstart+1:Pend-1),',');
        %Get the actual variables passed there
        for i=1:length(fcnInName)
            fcnInVar{i} = evalin('caller',fcnInName{i});
        end
    end
    %Check if the code in this line is relevant to validation
    useCode = iscoderelevant(String);
    %Check if the current line is defining one of the validation parameters
    for i=1:length(Parameters)
        SearchStr = [Parameters(i).name ' ='];
        SearchStr2 = [Parameters(i).name '='];
        if contains(String,SearchStr) || contains(String,SearchStr2)
            useCode = false;
        end
    end
    %If the code in the line is useful, then save it
    if useCode
        code{end+1} = String;
    end
end
if ~callerisfunction
    fcnInName = [];
    fcnInVar = [];
end

%-----------------------------------------------------
% Code-writter
%-----------------------------------------------------

%Create a temporal file to contain the script as a function
FID = fopen('process2validate.m','w');
%Start including the function definition header
fprintf(FID, ['function ',OutParam,' = process2validate(varargin)\n']);
%Unpack the varargin inputs
fprintf(FID,'varargin = varargin{1}; \n');
%If the caller was a function, make sure the original inputs are passed
if callerisfunction
    for j=1:length(fcnInName)
        fprintf(FID,[fcnInName{j} '= varargin{',num2str(j),'}; \n']);
    end
else
    j = 0;
end
%Ensure the validation parameters are extracted as well
for i=1:length(Parameters)
    fprintf(FID,[Parameters(i).name '= varargin{',num2str(i+j),'}; \n']);
end
%Insert the code from the user
for i=1:length(code)
    fprintf(FID,[code{i} '\n']);
end
%Close the file
fclose(FID);

%-----------------------------------------------------
% Run validation
%-----------------------------------------------------

%Check how many inputs will be passed
Noriginalinputs = length(fcnInName);
Nvarlinputs = size(validationParam,2);
%If the function requires stable inputs, set them now in the varargin
if callerisfunction
    argsin(1:Noriginalinputs) = fcnInVar;
else
    argsin = {};
end
%Prepare the output variable container
out = [];
%Run over all validation parameter permutations
for i=1:length(validationParam)
    argsin(Noriginalinputs+1:Noriginalinputs+Nvarlinputs) = validationParam(i,:);
    %Run the build function
    out(end+1,:) = process2validate(argsin);
    %Calculate status of validation statistics
    meanOut = mean(out,1);
    stdOut = std(out,[],1);
    %If user passes optional plotting hook, then prepare the plot
    if ~isempty(AxisHandle)
        cla(AxisHandle)
        Ax  = 1:length(meanOut);
        plot(AxisHandle,Ax,meanOut,'k','LineWidth',1)
        hold(AxisHandle,'on')
        f = fill(AxisHandle,[Ax fliplr(Ax)] ,[meanOut+stdOut fliplr(meanOut-stdOut)],...
            'b','LineStyle','none');
        f.FaceAlpha = 0.5;
        hold(AxisHandle,'off')
        axis(AxisHandle,'tight')
        grid(AxisHandle,'on')
        box(AxisHandle,'on')
        title(sprintf('Run %i/%i',i,length(validationParam)))
        drawnow
    end
    meanOut =meanOut.';
    stdOut = stdOut.';
end

%Delete the temporal file
delete('process2validate.m')


end









