function mfield_out = mask_field( x, mask, xvals_vecs, asnan )
% NEWFUN serves as a function template.
%--------------------------------------------------------------------------
% ARGUMENTS
% x
% mask
% xvals_vecs
% asnan
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if nargin < 4
    asnan = 1;
end

RK = @(y)boxker(y, 0.5, 1);
if nargin < 3
    mfield = @(tval) applyconvfield(tval, nan2zero(mask), RK, 1);
else
    mfield = @(tval) applyconvfield(tval, nan2zero(mask), RK, 1, xvals_vecs);
end

if ~asnan
    mfield_out = inf2zero(mfield(x)./mfield(x));
else
    mfield_out = inf2nan(mfield(x)./mfield(x)); %Returns NaN outside of the mask!
end

end

