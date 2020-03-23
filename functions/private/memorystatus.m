%
% MEMORYSTATUS Get memory information from OS
%
% [memfree,memsize] = MEMORYSTATUS()
%  Returns the available memory to the MATLAB process (memfree) in GB and
%  the total physical memory (RAM) available to the OS (memsize) in GB.
%  The function works for all OS systems: WinOS, UNIX and MacOS.
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function [memfree,memsize] = memorystatus()


if ismac
    
    % Get PID of MATLAB process
    [~,str] = system(['ps -o ppid,command |grep ',version('-release'),'.app']);
    slash = strfind('/',str);
    PID = verStr(1:slash(1)-1);
    % Get memory information of MATLAB process
    [~,info] = system(['top -pid ',PID,' -l 1 -s 0 | grep PhysMem']);
    % Parse results
    pos = strfind(info,' ');
    memsize = info(pos(1)+1:pos(2)-2);
    memsize = str2double(memsize)/1e3; %GB
    memfree = info(pos(end-1)+1:pos(end)-2);
    memfree = str2double(memfree)/1e3; %GB
    
elseif isunix
    
    % Retrieve free RAM info via OS "free" command
    [~,output] = unix('free | grep Mem');
    % Parse results
    stats = str2double(regexp(output, '[0-9]*', 'match'));
    memsize = stats(1)/1e6; %GB
    memfree = (stats(3)+stats(end))/1e6; %GB
    
elseif ispc
    
    % Retrieve system memory info via built-in function
    [~,sys] = memory;
    memsize = sys.PhysicalMemory.Total/1e9; %GB
    memfree = min(sys.PhysicalMemory.Available,sys.VirtualAddressSpace.Available)/1e9; %GB
    
else
    disp('OS platform not supported')
end

end