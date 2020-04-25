    % Tests whether functions from the Optimization Toolbox can be used

function tf = optimtoolbox_installed()

licensed = license('test','optimization_toolbox');
installed = exist('lsqnonlin','file') && exist('optimoptions','file');

tf = licensed && installed;

return
