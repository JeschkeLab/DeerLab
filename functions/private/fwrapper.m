function [f,J,ctrlout,evalJ]=fwrapper(x,ctrlin,func,evalJ,Joption,varargin)
% wrapper for elsunc compatibility with Matlab
%evalJ=varargin{1};Joption=varargin{2};
x=x(:);
if ctrlin==1||ctrlin==-1% computes f
    J=[];
    try
        f=func(x,varargin{:});ctrlout=ctrlin;
    catch% f is not computable at x
        if ctrlin==1
            ctrlout=-1;
        else
            ctrlout=-20;
        end
        f=[];
    end
elseif ctrlin==2&&Joption==1% computes analytical J
    try
        [f,J]=func(x,varargin{:});ctrlout=ctrlin;evalJ=evalJ+1;
    catch
        ctrlout=-30;% J is not computable at x
        f=[];J=[];
    end
else
    ctrlout=0;% J wil be computed by differences
    f=[];J=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    