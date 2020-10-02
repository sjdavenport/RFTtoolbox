clear all
close all

% Changes here for playing around
D   = 3;
pad = 0;
resadd = [ 0, 1, 2, 3 ];
enlarge = zeros( [ 1 length( resadd ) ] );
%enlarge = ceil( resadd / 2 );

% get dimension of non padded domain
switch D
    case 1
        dim  = [ 20, 1 ];
        dimp = [ 20 + 2 * pad, 1 ];
    case 2
        dim  = [ 3, 4 ];
        dimp = dim + 2 * pad * [ 1, 1 ];
    case 3
        dim = [ 4, 4, 4 ];
        dimp = dim + 2 * pad * [ 1, 1, 1 ];        
end

xvals  = cell( [1 D] );
for d = 1:D
    xvals{d} = 1:dimp(d); 
end


% Mask

    mask = true( dim );
    mask = pad_vals( mask, pad, false );
    mask( 1, 1 ) = false;
    %mask = logical( pad_vals( mask, pad) );
    mask(2,1) = false;

vol = NaN * ones( length( resadd ), D );

% voxel manifolds for two resolutions
for r = 1:length( resadd )
    % high resolution mask
    mask_hr      = mask_highres( mask, resadd(r), enlarge(r) );
    % voxel manifold from high resolution mask with euclidean metric
    voxmfd       = VoxManifold( mask_hr );
    % correct the xvals vector to high resolution
    voxmfd.xvals = xvals_highres( xvals, resadd(r), enlarge(r) );
    vol(r,:)     = LKC_est( voxmfd );
end

% Volume of the domain and of the boundary
vol

