mask = ones(10,10);
bdry_voxels( logical(mask), 'full' )
bdry_voxels_mod( logical(mask), 'full' )

%%
bdry_voxels_mod( logical(mask), 'x' )
bdry_voxels( logical(mask), 'x' )

