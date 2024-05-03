function weights = getweights( mask_hr )
% GETWEIGHTS( mask_hr )
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% 	mask_hr     the high resolution mask
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% mask_hr = ones(3); mask_hr(3,3) = 0;
% getweights( mask_hr )
% %% 2D Weights
% % resolution added
% resAdd  = 1;
% % create a mask and show it
% Sig  = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
% mask = logical( Sig > 0.02 & Sig < 1.1 );
% 
% % enlarged domain
% [mask_hr, weights, old_weights] = mask_highres( mask, resAdd, 1 );
% subplot(2,1,1); imagesc(weights); title('New Weights')
% subplot(2,1,2); imagesc(old_weights); title('Old Weights')
% 
% %% 3D weights (the examples we considered today :))
% example_mask = ones(3,3,3);
% example_mask(1,1,1) = 0;
% 
% weights = getweights(example_mask)
% weights(2,2,2) % I.e. 7/8 like I said haha
% 
% % Your example (I think this was your example)
% example_mask = ones(3,3,3);
% example_mask(1,1,2) = 0;
% 
% weights = getweights(example_mask)
% weights(2,2,2) 
%
% %% Testing that the volume works correctly!
% get_volume = @( weights, resadd, D ) sum( weights(:) ) * 1/(resadd+1)^D;
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport, Fabian Telschow
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

% Get size of the mask
Dim_hr = size(mask_hr);

% Get Dimension of the mask
D = length(Dim_hr);
if D == 2 && ( Dim_hr(1) == 1 || Dim_hr(2) == 1 )
    D = 1;
end

%%  Main Function Loop
%--------------------------------------------------------------------------

% Easy fix for the bug. Shown in test_volume_highres_mask.
% Maybe there is something more smart possible. To see the bug uncomment
% line 65:69. Also uncomment line 98/99

% Fill mask into mask enlarged by padding 0's
[ mask_hr, locs ] = pad_vals( mask_hr );

% Get new resolution
Dim_hr = size(mask_hr);

% Divide voxels in the high res mask in two in each dimension 
mask_with_divided_voxels = mask_highres_divide( mask_hr, 2 );

% Erode the new mask by one voxel
dilatedanddividedmask = dilate_mask( mask_with_divided_voxels, -1 );
clear mask_with_divided_voxels

% Sum up the number of small voxels within each larger voxel
ones_array = ones( ones( 1, D ) * 2 );
sum_within_each_voxel_large = convn( dilatedanddividedmask, ones_array );
clear dilatedanddividedmask

% The above has too much information due to the convolution we only want
% the corner sum. Then we need to divide by the total number of small
% voxels within each large voxel which is 2^D.
index = cell( [ 1 D ] );
for d = 1:D
    index{d} = 2:2:(2*Dim_hr(d));
end
weights = sum_within_each_voxel_large( index{:} ) / 2^D;

% Remove the extra padded 0 values
weights = weights(locs{:});
end

