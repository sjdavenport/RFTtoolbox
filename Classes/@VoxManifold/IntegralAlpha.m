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
%            pi/2 or 3*pi/2 (or pi)
%--------------------------------------------------------------------------
% OUTPUT
%   angles   the inside angles between the voxel manifold and the orthogonal
%            plane to the edge (wrt voxmfd.g)  
%--------------------------------------------------------------------------
% DEVELOPER TODOs:
% This code might be slower than possible as many zeros are carried around.
% We should consider at some point to speed it up, but for now it at least
% provides the correct results!
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
            padz  = [1 0 0];
        case 'y'
            direc = [ 2 3 1 ];
            padz  = [0 1 0];
        case 'z'
            direc = [ 3 1 2 ];
            padz  = [0 0 1];
    end

    % Create an orthonormal frame containing the chosen direction as one vector
    [ ~, b2, b3 ] = OrthNormFrame( voxmfd, direc, mask.(type), 0 );
    
    % Normal of plane defined by b2 and b3
    n = cross( b2, b3 );
    n = n ./ sqrt( n.' * n );

    n1 = n(:,:,:,1);
    n2 = n(:,:,:,2);
    n3 = n(:,:,:,3); 

    % Get the cutting angle between the plane spanned by b2 and b3 and the
    % voxelmanifold
    switch type
        case 'x'
            angle = pi - acos( n2 .* n3 ./ sqrt( n1.^2 + n2.^2 ) ...
                                   ./ sqrt( n1.^2 + n3.^2 ) );
        case 'y'
            angle = pi - acos( n1 .* n3 ./ sqrt( n1.^2 + n2.^2 ) ...
                                   ./ sqrt( n2.^2 + n3.^2 ) );
        case 'z'
            angle = pi - acos( n2 .* n1 ./ sqrt( n1.^2 + n3.^2 ) ...
                                   ./ sqrt( n2.^2 + n3.^2 ) );
    end

    % NaN to 0
    angle = angle.field;
    angle(isnan(angle)) = 0;
    angle = pad_vals(angle, padz);

    % pad a 0 to the angle field
    NN = 1:length(angle(:));

    eucangle.(type) = pad_vals(eucangle.(type), padz);

    if type == 'y'
        eucangle.(type) = permute(eucangle.(type), [2 1 3] );
        angle = permute(angle, [2,1,3] );
    elseif type == 'z'
        eucangle.(type) = permute(eucangle.(type), [3 2 1] );
        angle = permute(angle, [3 2 1] );
    end

    % convex edge
    angle( eucangle.(type) == 7 ) = -2*(pi - angle( eucangle.(type) == 7 ));% 0; %angle( eucangle.(type) == 7 ) * 2;%
    % convex edge meeting double convex edge
    angle( eucangle.(type) == 8 )  = -angle( eucangle.(type) == 8 );
    % convex edge meeting concave edge
    angle( eucangle.(type) == 42 ) = (angle( NN(eucangle.(type) == 42) - 1 )...
                                      + angle( NN(eucangle.(type) == 42) + 1 )) / 2;

    angles.(type) = angle(2:(end-1), :, :);

    if type == 'y'
         angles.(type) = permute(angles.(type), [2 1 3] );
    elseif type == 'z'
        angles.(type) = permute(angles.(type), [3 2 1] );
    end
end

return