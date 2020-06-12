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
% 
% %% View masked field (old version)
% xvals_fine = 0.5:0.1:3.5;
% xvaluesatvoxels = xvals2voxels(xvals_fine,2)
% masked_field = @(x) mask_field(x, mask).*field(x);
% fieldeval = reshape(masked_field(xvaluesatvoxels), [length(xvals_fine),length(xvals_fine)]);
% surf(fieldeval)
% 
% Y = [1,1,1;10,1,1;1,1,1] %I.e so the peak will be outside the mask!
% mask = [0,1,1;0,1,1;0,1,1];
% findconvpeaks(Y, 3, [3,3]', mask) %Need to fix so that initialization at [2,2] is fine as well!
% 
% %% Display masked field using applyconvfield
% field = @(tval) applyconvfield(tval, Y, FWHM, -1, xvals, mask);
% masked_field = @(x) mask_field(x, mask).*field(x);
% xvals_fine = 0.5:0.1:(length(Y)+0.5);
% xvaluesatvoxels = xvals2voxels(xvals_fine,2);
% fieldeval = reshape(masked_field(xvaluesatvoxels), [length(xvals_fine),length(xvals_fine)]);
% surf(fieldeval)
% 
% %% View masked field (old version)
% FWHM = 3; D = 2;
% Y = [10,1,1;1,1,1;10,1,1]; %I.e. so the peak will be outside the mask!
% mask = logical([1,1,1;0,1,1;1,1,1]);
% xvals_fine = 0.5:0.1:3.5;
% xvaluesatvoxels = xvals2voxels(xvals_fine,2)
% field = @(tval) applyconvfield(tval, Y, FWHM, -1, 1:3, mask);
% 
% masked_field = @(x) mask_field(x, mask).*field(x);
% fieldeval = reshape(masked_field(xvaluesatvoxels), [length(xvals_fine),length(xvals_fine)]);
% surf(fieldeval)
% 
% Y = [1,1,1;10,1,1;1,1,1] %I.e so the peak will be outside the mask!
% mask = [0,1,1;0,1,1;0,1,1];
% findconvpeaks(Y, 3, [2,2]', mask) %Need to fix so that initialization at [2,2] is fine as well!
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

