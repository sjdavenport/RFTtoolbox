function bdry_params = getBdryparams( field, c )
% getBdryparams( field, c ) finds edges and weights needed to linearly
% interpolate a discretized field to find the locations where the true
% continuous field = c
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% field     an array representing a scalar field over a domain.
% c:        threshold value for excursion
%--------------------------------------------------------------------------
% OUTPUT
%   bdry_params  a struct containing the edge locations either side of the
%                true continuous boundary where the field = c, as well as
%                the weights that the adjacent voxels need to be multiplied
%                by to give the precise location where the field = c
%                assuming that the gradient of the signal is linear between
%                adjacent voxels. 
%
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR:  Fabian Telschow, Alex Bowring
%--------------------------------------------------------------------------
%% Check mandatory input and get important constants
%--------------------------------------------------------------------------


%% Add/check optional values
%--------------------------------------------------------------------------

A_c = ( field >= c );
dim = size( A_c );
D   = length( dim );

switch D
    case 2
        %% Case 2D random field with 4-connectivity
        %----------------------------------------------------------------------
        %%% Horizontal edges (note that this uses the matlab image nomenclature,
        %%% usually the second component might be called vertical)
        horz = A_c( :, 2:end ) | A_c( :, 1:end-1 );
        % Compute the left shifted horizontal edges
        lshift            = A_c; % initialize
        lshift(:,1:end-1) = horz;
        lshift            = lshift & ~A_c;
        
        % Compute the right shifted horizontal edges
        rshift          = A_c; % initialize
        rshift(:,2:end) = horz;
        rshift          = rshift & ~A_c;
        
        %%% Vertical edges (note that this uses the matlab image nomenclature,
        %%% usually the first component might be called horizontal)
        vert = A_c( 1:end-1, : ) | A_c( 2:end, : );
        % Compute the right shifted horizontal edges
        ushift            = A_c;
        ushift(1:end-1,:) = vert;
        ushift            = ushift & ~A_c;
        
        % Compute the down shifted vertical edges
        dshift          = A_c;
        dshift(2:end,:) = vert;
        dshift          = dshift & ~A_c;
        
        % Computing the weights for the weighted linear boundary of the
        % field
        lshift_w1 = abs( field( lshift( :, [ dim(2) 1:dim(2)-1 ] ) ) - c )...
                        ./abs( field( lshift ) - field( lshift( :, [ dim(2) 1:dim(2)-1 ] ) ) );
        lshift_w2 = abs( field( lshift) - c)...
                        ./abs( field( lshift ) - field( lshift( :, [ dim(2) 1:dim(2)-1 ] ) ) );

        rshift_w1 = abs( field( rshift( :, [ 2:dim(2) 1 ] ) ) - c )...
                        ./ abs( field( rshift ) - field( rshift( :, [ 2:dim(2) 1 ] ) ) );
        rshift_w2 = abs( field( rshift) - c)...
                        ./ abs( field( rshift ) - field( rshift( :, [ 2:dim(2) 1 ] ) ) );

        ushift_w1 = abs( field( ushift( [ dim(1) 1:dim(1)-1 ], : ) ) - c )...
                        ./ abs( field( ushift ) - field( ushift( [ dim(1) 1:dim(1)-1 ], : ) ) );
        ushift_w2 = abs( field( ushift ) - c )...
                        ./ abs( field( ushift ) - field( ushift( [ dim(1) 1:dim(1)-1 ], : ) ) );

        dshift_w1 = abs( field( dshift( [ 2:dim(1) 1 ], : ) ) - c )...
                        ./ abs( field( dshift ) - field( dshift( [ 2:dim(1) 1 ], : ) ) );
        dshift_w2 = abs( field( dshift ) - c )...
                        ./ abs( field( dshift ) - field( dshift( [ 2:dim(1) 1 ], : ) ) );
        
        % Compute the length of the boundary
        len = length( lshift_w1 ) + length( rshift_w1 ) +...
              length( ushift_w1 ) + length( dshift_w1 );
        
        % Storing parameters in a structure
        bdry_params = struct( 'length', len, ...
                              'lshift', struct( 'edges', lshift,...
                                                'w1', lshift_w1,...
                                                'w2', lshift_w2 ), ...
                              'rshift', struct( 'edges', rshift,...
                                                'w1', rshift_w1,...
                                                'w2', rshift_w2 ), ...
                              'ushift', struct( 'edges', ushift,...
                                                'w1', ushift_w1,...
                                                'w2', ushift_w2 ), ...
                              'dshift', struct( 'edges', dshift,...
                                                'w1', dshift_w1,...
                                                'w2', dshift_w2 ) );
                         
    case 3
        %% Case 3D random field with 6-connectivity
        %------------------------------------------------------------------                      
        %%% Horizontal edges (note that this uses the matlab image nomenclature, 
        %%% usually the second component might be called vertical)
        horz = A_c( :, 2:end, : ) | A_c( :, 1:end-1, : );
        % Compute the left shifted horizontal edges
        lshift              = A_c; % initialize
        lshift(:,1:end-1,:) = horz;
        lshift              = lshift & ~A_c;
        % Compute the right shifted horizontal edges
        rshift            = A_c; % initialize
        rshift(:,2:end,:) = horz;
        rshift            = rshift & ~A_c;
        
        %%% Vertical edges (note that this uses the matlab image nomenclature,
        %%% usually the first component might be called horizontal)
        vert = A_c( 1:end-1, :, : ) | A_c( 2:end, :, : );
        % Compute the right shifted horizontal edges
        ushift              = A_c;
        ushift(1:end-1,:,:) = vert;
        ushift              = ushift & ~A_c;
        % Compute the down shifted vertical edges
        dshift            = A_c;
        dshift(2:end,:,:) = vert;
        dshift            = dshift & ~A_c;
        
        %%% depth edges
        depth = A_c( :, :, 1:end-1 ) | A_c( :, :, 2:end );
        % Compute the back shifted depth edges
        bshift               = A_c;
        bshift(:,:,1:end-1 ) = depth;
        bshift               = bshift & ~A_c;
        % Compute the fron shifted depth edges
        fshift              = A_c;
        fshift(:,:,2:end)   = depth;
        fshift              = fshift & ~A_c;

        % Computing the weights for the weighted linear boundary of the
        % field
        lshift_w1 = abs( field( lshift( :, [ dim(2) 1:dim(2)-1 ], : ) ) - c )...
                        ./ abs( field( lshift ) - field( lshift( :, [ dim(2) 1:dim(2)-1 ], : ) ) );
        lshift_w2 = abs( field( lshift ) - c )...
                        ./ abs( field( lshift ) - field( lshift( :, [ dim(2) 1:dim(2)-1 ], : ) ) );

        rshift_w1 = abs( field( rshift( :, [ 2:dim(2) 1 ], : ) ) - c )...
                        ./ abs( field( rshift ) - field( rshift( :, [ 2:dim(2) 1 ], : ) ) );
        rshift_w2 = abs( field( rshift ) - c)...
                        ./ abs( field( rshift ) - field( rshift( :, [ 2:dim(2) 1 ], : ) ) );

        ushift_w1 = abs( field( ushift( [ dim(1) 1:dim(1)-1 ], :, : ) ) - c )...
                        ./ abs( field( ushift ) - field( ushift( [ dim(1) 1:dim(1)-1 ], :, : ) ) );
        ushift_w2 = abs( field( ushift ) - c)...
                        ./ abs( field( ushift ) - field( ushift( [ dim(1) 1:dim(1)-1 ], :, : ) ) );

        dshift_w1 = abs( field( dshift( [ 2:dim(1) 1 ], :, : ) ) - c )...
                        ./ abs( field( dshift ) - field( dshift( [ 2:dim(1) 1 ], :, : ) ) );
        dshift_w2 = abs( field( dshift ) - c)...
                        ./ abs( field( dshift ) - field( dshift( [ 2:dim(1) 1 ], :, : ) ) );

        bshift_w1 = abs( field( bshift( :, :, [ dim(3) 1:dim(3)-1 ] ) ) - c )...
                        ./ abs( field( bshift ) - field( bshift( :, :, [ dim(3) 1:dim(3)-1 ] ) ) );
        bshift_w2 = abs( field( bshift ) - c)...
                        ./abs( field( bshift ) - field( bshift( :, :, [ dim(3) 1:dim(3)-1 ] ) ) );

        fshift_w1 = abs( field( fshift( :, :, [ 2:dim(3) 1 ] ) ) - c )...
                        ./ abs( field( fshift ) - field( fshift( :, :, [ 2:dim(3) 1 ] ) ) );
        fshift_w2 = abs( field( fshift ) - c )...
                        ./ abs( field( fshift ) - field( fshift( :, :, [ 2:dim(3) 1 ] ) ) );
        
        % Compute the length of the boundary
        len = length( lshift_w1 ) + length( rshift_w1 ) + ...
              length( ushift_w1 ) + length( dshift_w1 ) + ...
              length( bshift_w1 ) + length( fshift_w1 );

        
        % Create structure for storing parameters
        bdry_params = struct( 'length', len, ...
                              'lshift', struct( 'edges', lshift,...
                                                'w1', lshift_w1,...
                                                'w2', lshift_w2 ), ...
                              'rshift', struct( 'edges', rshift,...
                                                'w1', rshift_w1,...
                                                'w2', rshift_w2 ), ...
                              'ushift', struct( 'edges', ushift,...
                                                'w1', ushift_w1,...
                                                'w2', ushift_w2 ), ...
                              'dshift', struct( 'edges', dshift ,...
                                                'w1', dshift_w1,...
                                                'w2', dshift_w2 ), ...
                              'bshift', struct( 'edges', bshift ,...
                                                'w1', bshift_w1,...
                                                'w2', bshift_w2 ), ...
                              'fshift', struct( 'edges', fshift,...
                                                'w1', fshift_w1,...
                                                'w2', fshift_w2 ) );        
end