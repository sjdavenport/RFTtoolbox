function [mfield_out, ceq] = mask_field( x, mask, asnan, xvals_vecs )
% MASK_FIELD( x, mask, xvals_vecs, asnan ) generates an indicator of
% whether a given point x is within 0.5 distance of a voxel that lies in
% the mask. Note that mask_field (like the applyconvfields) is designed to
% be used locally. If you would like to generate a high resolution mask 
% over the whole field you should use mask_highres
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  x           a D by nvals vector 
%  mask        a D dimension lattice of 0/1 indicating where the mask is
% Optional
%  xvals_vecs    a D-dimensional cell array whose entries are vectors giving 
%               the xvalues at each each dimension along the lattice. It 
%               assumes a regular, rectangular lattice (though within a given
%               dimension the voxels can be spaced irregularly).
%               I.e suppose that your initial lattice grid is a
%               4by5 2D grid with 4 voxels in the x direction and 5 in
%               the y direction. And that the x-values take the values:
%               [1,2,3,4] and the y-values take the values: [0,2,4,6,8].
%               Then you would take xvals_vecs = {[1,2,3,4], [0,2,4,6,8]}.
%               The default is to assume that the spacing between the
%               voxels is 1. If only one xval_vec direction is set the
%               others are taken to range up from 1 with increment given by
%               the set direction.
%  asnan        0/1 whether to returns -Infs or NaNs outside of the mask.
%               Default is 1 i.e. to return NaNs.
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if ~exist('asnan', 'var')
    asnan = 1;
end

% Let RK be the box kernel
RK = @(y)boxker(y, 0.5, 1);

% Calulate the masked_field handle (if xvals_vecs is not an input don't use
% it when computing the field)
if ~exist('xvals_vecs', 'var')
    mfield = @(tval) applyconvfield(tval, nan2zero(mask), RK, NaN, 1);
else
    mfield = @(tval) applyconvfield(tval, nan2zero(mask), RK, NaN, 1, ...
                                                               xvals_vecs);
end

% Divide to normalize (as at edges of the box kernel voxels can be counted
% twice)
val = mfield(x);
if asnan == -1
    mfield_out = inf2nan(val./val);
    mfield_out(isnan(mfield_out)) = -1;
    mfield_out = -mfield_out;
else
    if isinf(asnan)
%         mfield_out = mfield(x)./mfield(x);
%         mfield_out(isnan(mfield_out)) = 0.2;
        mfield_out = nan2inf(val./val,-1);
    elseif asnan == 0
        mfield_out = val./val;
        mfield_out = nan2zero(mfield_out);
        mfield_out = inf2zero(mfield_out);
    elseif asnan
        mfield_out = inf2nan(val./val); %Returns NaN outside of the mask!
    end
end

ceq = [];

end

