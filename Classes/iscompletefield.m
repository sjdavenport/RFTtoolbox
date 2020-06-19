function out = iscompletefield( var )
% iscompletefield( var ) tests whether the input is of Class Field and is a
% complete field, i.e. field, mask and xvals are properly specified.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  var   any varible type
%
%--------------------------------------------------------------------------
% OUTPUT
% out  a logical which is true, if 'var' is of Class Field and complete.
%      false otherwise.
%
%--------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% Author: Fabian Telschow
%--------------------------------------------------------------------------


%% Main function
%--------------------------------------------------------------------------
out = 0;

if isa( var, 'Field' )
    if var.complete
        out = true;    
    end
end

return