function [ b1, b2, b3 ] = OrthNormFrame( voxmfd, direc, mask, vectorized )
% Onf( voxmfd, direc, masked ) This function computes boundary voxels and weights for a given mask. Note
% that the weights only make sense for masks which are resolution increased
% by an odd number using mask_highres.m.
% The weights are used in LKC estimation and assume trapozoidal
% integration.
% The boundary is splitted into several subparts representing the faces and
% edges for specific fixed coordinates.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   voxmfd     a logical T_1 x ... x T_D array.
% Optional
%   direc   a 1 x D numeric array containing the indices of the standard
%           euclidean basis and how it is used.
%   mask    a logical array
%   vectorized  a logical, if true the output fields domain is vectorized
%               for memory efficiency
%--------------------------------------------------------------------------
% OUTPUT
%   ONF     an object of type field containing an orthormal frame with
%           respect to g in the fiber. 
%--------------------------------------------------------------------------
% DEVELOPER TODOS: 
%--------------------------------------------------------------------------
% EXAMPLES
%--------------------------------------------------------------------------
% AUTHORS: Fabian Telschow
%--------------------------------------------------------------------------

%% Check input and get important constants from the mandatory input
%--------------------------------------------------------------------------

% Get dimension of the voxel manifold
D = voxmfd.D;

% Reject input if D > 3
if D > 3
    error( "The method is currently only implemented for D < 4." )
end

%%  Add/check optional values
%--------------------------------------------------------------------------

if ~exist( 'direc', 'var' )
    direc = [ 1 2 ];
    if D == 3
        direc = [ 1 2 3 ];
    end
end

if ~exist( 'vectorized', 'var' )
    vectorized = 0;
end

%% Main function
%--------------------------------------------------------------------------
        
switch D
    case 2
    case 3
        
        % Get the scalar field entries from the Riemannian metric
        gkk = Subfield( voxmfd.g, { ':', ':', ':', direc(1), direc(1) } );
        gkl = Subfield( voxmfd.g, { ':', ':', ':', direc(1), direc(2) } );
        gll = Subfield( voxmfd.g, { ':', ':', ':', direc(2), direc(2) } );
        
        if ~all( mask(:) )
            gkk = Subfield( collapse( gkk ), { mask(:) } );
            gkl = Subfield( collapse( gkl ), { mask(:) } ); 
            gll = Subfield( collapse( gll ), { mask(:) } );
            g   = Subfield( collapse( voxmfd.g ), { mask(:), ':', ':' } );
        end

            % Get the constant euclidean basis fields
            Ek = constfield( sbasis( direc(1), D ), gll.mask, gll.xvals );
            El = constfield( sbasis( direc(2), D ), gll.mask, gll.xvals );
            Es = constfield( sbasis( direc(3), D ), gll.mask, gll.xvals );

            %%% Get the basis vectors
            % Basis field along first direction b1 = e_k / g_kk
            b1 = gkk.^(-1/2) .* Ek;

            % compute a normalization constant
            c = sqrt( gkk .* gll - gkl.^2 );

            % Basis field along second direction b2 = e_k / g_kk
            b2 = gkk.^(1/2) .* El - gkk.^(-1/2) .* gkl .* Ek;
            b2 = b2 ./ c;

            % Basis field along third direction
            % b3 = G^-1 e_s / sqrt( e_s^T G^-1 e_s )
            b3 = invsym( g ) * squeeze( Es );
            b3 = b3 ./ sqrt( Es.' * b3 );
            
            % Reshape basis vectors to have correct size
            if ~all( mask(:) ) && ~vectorized
                tmpfiber   = zeros( [ prod( voxmfd.size ), 3 ] );
                tmpfiber( mask(:), : ) = b1.field;
                b1 = voxmfd.g;
                b1.mask = mask;
                b1.field = reshape( tmpfiber, [ voxmfd.size, 3 ] );
                
                tmpfiber   = zeros( [ prod( voxmfd.size ), 3 ] );
                tmpfiber( mask(:), : ) = b2.field;
                b2 = voxmfd.g;
                b2.mask = mask;
                b2.field = reshape( tmpfiber, [ voxmfd.size, 3 ] );
                
                tmpfiber   = zeros( [ prod( voxmfd.size ), 3 ] );
                tmpfiber( mask(:), : ) = b3.field;
                b3 = voxmfd.g;
                b3.mask = mask;
                b3.field = reshape( tmpfiber, [ voxmfd.size, 3 ] );
            end
            
end

return