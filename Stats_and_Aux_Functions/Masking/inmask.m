function out = inmask(point, mask)
% INMASK(point, mask) tests whether a given point is in a given mask.
%--------------------------------------------------------------------------
% ARGUMENTS
% point: a D x 1 vector where D is the number of dimensions, representing 
%        the point being tested.
% mask: a logical array of size N_1 x N_2 x ... x N_D, where 
%       N_1, N_2, ..., N_D are the dimensions of the mask.
%--------------------------------------------------------------------------
% OUTPUT
% out: a logical scalar, true if the point is in the mask, false otherwise.
%--------------------------------------------------------------------------
% EXAMPLES
% %1D
% inmask( 0.9, ones(1,5))
% inmask( 0.45, ones(1,5))
% inmask( 5.45, ones(1,5))
% inmask( 5.55, [ones(1,5),0])
% 
% %2D
% mask = [0,0,1;0,0,1;0,0,1]
% inmask([1,2]', mask)
% inmask([2.5,2.6]', mask)
% inmask([2.5,2.4]', mask)
%
% % 3D MNImask
% MNImask = imgload('MNImask')
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
[~,nearest_ones] = nearest_voxels( point );

s_mask = size(mask);

Dim = s_mask(s_mask>1);

n_nearestvoxels = size(nearest_ones,2);
nearestvoxels_inmask = zeros(1,n_nearestvoxels);

for I = 1:size(nearest_ones,2)
    if any(nearest_ones(:,I)<= 0) || any(nearest_ones(:,I)> Dim')
        nearestvoxels_inmask(I) = 0;
    else
        idx = convind(nearest_ones(:,I),s_mask);
        nearestvoxels_inmask(I) = mask(idx);
    end
end

if any(nearestvoxels_inmask)
    out = 1;
else
    out = 0;
end

end
