function [ lowerbounds, upperbounds ] = assign_bounds( points, mfield )
% ASSIGN_BOUNDS( points, bfield )
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  points
%  bfield
%--------------------------------------------------------------------------
% OUTPUT
%
%--------------------------------------------------------------------------
% EXAMPLES
% boundary = ones(3,3); boundary(2,2) = 0;
% bfield =  @(x) mask_field( x, boundary, 0 );
% [lb,ub] = assign_bounds( [2,2.5]', bfield )
% [lb,ub] = assign_bounds( [3,3]', bfield )
% [lb,ub] = assign_bounds( [0.5,0.5]', bfield )
%
% boundary = ones(2); boundary(1,1) = 0;
% bfield =  @(x) mask_field( x, boundary, 0 );
% [lb,ub] = assign_bounds( [1.5,1.5]', bfield )
%
% boundary = ones([2,2,2]); boundary(1,1,1) = 0;
% bfield =  @(x) mask_field( x, boundary, 0 );
% [lb,ub] = assign_bounds( [1.5,1.5,1.5]', bfield )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'opt1', 'var' )
    % default option of opt1
    opt1 = 0;
end

% if isnumeric(bfield) || islogical(bfield)
%     bfield =  @(x) mask_field( x, boundary );
% end

%%  Main Function Loop
%--------------------------------------------------------------------------
npeaks = size(points,2);

lowerbounds = cell(1,npeaks);
upperbounds = cell(1,npeaks);
for I = 1:npeaks
    point = points(:,I);
    [ nearby_voxels, nearest_ones ] = nearest_voxels( point );
    nearby_voxels_inmask_ind = mfield(nearby_voxels);
    nearest_ones_inmask_ind = mfield(nearest_ones);
%     totalinmask = sum(nearby_voxels_inmask_ind);
    
    if all(nearby_voxels_inmask_ind)
        lowerbounds{I} = floor(point) - 0.5;
        upperbounds{I} = ceil(point) + 0.5;
    elseif ~any(nearest_ones_inmask_ind)
        error('The original points must lie in the mask!')
    else
        voxelsinmask = nearest_ones(:, logical(nearest_ones_inmask_ind));
        lowerbounds{I} = voxelsinmask - 0.5;
        upperbounds{I} = voxelsinmask + 0.5;
    end
end

% By searching for voxels that share all values except in one direction can
% sparsify the lower and upper bounds!
end

