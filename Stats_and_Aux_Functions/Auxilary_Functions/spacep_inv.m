function point = spacep_inv(point_out, resadd)
% INVERSE_SPACEP(POINT_OUT, RESADD) converts the indices of a point in a 
% convolution field back to its original location using the given resadd value.
%--------------------------------------------------------------------------
% ARGUMENTS
% point_out  the indices of the point
% resadd     the resadd value used in the original call to spacep
%--------------------------------------------------------------------------
% OUTPUT
% point      the original location of the point
%--------------------------------------------------------------------------
% EXAMPLES
% enlarged_version = spacep( [3,3,3]', 20 )
% spacep_inv(enlarged_version, 20)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

spacing = 1/(1+resadd);
inverse_spacing = floor(1/spacing);
mod_spacing = 1/inverse_spacing;
point = (point_out - 1) * mod_spacing + 1;

end
