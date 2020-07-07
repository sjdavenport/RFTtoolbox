function angles = IntegralAlpha( voxmfd, mask, eucangle )
% IntegralAlpha( voxmfd, mask, type ) computes a field .
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  voxmfd    an object of class VoxManifold
%  mask
%  eucangle  the euclidean inside opening angles at the edges,
%            since a voxel manifold is built out of cubes these are either
%            pi/2 or 3*pi/2
%--------------------------------------------------------------------------
% OUTPUT
%   angles   the inside angles between the voxel manifold and the orthogonal
%            plane to the edge (wrt voxmfd.g)  
%--------------------------------------------------------------------------
% DEVELOPER TODOs:
%--------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow
%--------------------------------------------------------------------------

%% Main function
%--------------------------------------------------------------------------
angles = struct();

for type = ["x", "y", "z"]
    switch type
        case 'x'
            direc = [ 1 2 3 ];
        case 'y'
            direc = [ 2 3 1 ];
        case 'z'
            direc = [ 3 1 2 ];
    end

    % Create an orthonormal frame containing the chosen direction as one vector
    [ ~, b2, b3 ] = OrthNormFrame( voxmfd, direc, mask.(type), 1 );

    % Get the euclidean opening angles greater than pi/2
    greater_pi = eucangle.(type)( mask.(type)(:) ) > pi;

    % Normal of plane defined by b2 and b3
    n = cross( b2, b3 );
    n = n ./ sqrt( n.' * n );

    %
    n1 = Subfield( n, {':', 1} );
    n2 = Subfield( n, {':', 2} );
    n3 = Subfield( n, {':', 3} );

    % Get the cutting angle between the plane spanned by b2 and b3 and the
    % voxelmanifold
    switch type
        case 'x'
            angle = acos( n2 .* n3 ./ sqrt( n1.^2 + n2.^2 ) ...
                                   ./ sqrt( n1.^2 + n3.^2 ) );
        case 'y'
            angle = acos( n1 .* n3 ./ sqrt( n1.^2 + n2.^2 ) ...
                                   ./ sqrt( n2.^2 + n3.^2 ) );
        case 'z'
            angle = acos( n2 .* n1 ./ sqrt( n1.^2 + n3.^2 ) ...
                                   ./ sqrt( n2.^2 + n3.^2 ) );
    end

    % adjust for large angles
    angle.field( greater_pi ) = -angle.field( greater_pi );

    angles.(type) = angle.field;
end

return