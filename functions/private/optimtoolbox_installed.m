% Tests whether functions from the Optimization Toolbox can be used

function tf = optimtoolbox_installed()

licensed = license('test','optimization_toolbox');
installed = ~isempty(which('lsqnonlin')) && ~isempty(which('optimoptions'));

tf = licensed & installed;

return
