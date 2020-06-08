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
% FWHM = 3; 
% xvals_fine = 1:0.1:3; xvals_vecs = 1:3;
% xvaluesatvoxels = xvals2voxels({xvals_fine,xvals_fine});
% 
% Y = [10,1,1;1,1,1;10,1,1]; %I.e. so the peak will be outside the mask!
% mask = [1,1,1;0,1,1;1,1,1];
% 
% field = @(tval) applyconvfield(tval, Y, FWHM, -1, xvals_vecs, mask);
% mfield = @(x) mask_field(x, mask);
% masked_field = @(x) mfield(x).*field(x);
% fieldeval = reshape(masked_field(xvaluesatvoxels), [length(xvals_fine),length(xvals_fine)]);
% surf(fieldeval)
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

