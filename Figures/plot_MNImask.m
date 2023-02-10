MNImask = imgload('MNImask');
% [ bounds, bounded_mask ] = mask_bounds( MNImask );
MNImask_boundary = bndry_voxels(logical(MNImask), "full");

%%
f = Field();
f.mask = logical(MNImask_boundary);
f.xvals = {1:91,1:109,1:91};
plot_voxmf( f )