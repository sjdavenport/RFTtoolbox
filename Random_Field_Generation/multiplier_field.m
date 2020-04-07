function data = multiplier_field( base_fields, nSubj )
% multiplier_field( base_fields, nSubj ) generates samples from a mean zero
% unit variance field given by the base_fields using Gaussian multipliers.
%--------------------------------------------------------------------------
% ARGUMENTS
% base_fields   an array of dimension Dim x nBase containing
%               the base_fields. nBase is the amount of base_fields.
% nSubj         The number of subjects.
% shape_of_array    0/1/2, default is 0. Determines the shape of the array.
%                   If 0 data is Dim by nSubj,
%                   If 1 data is nSubj by Dim,
%                   If 2 data is prod(Dim) by nSubj instead of
%                   nSubj by Dim. 
%                   If 3 data is nSubj by prod(Dim) instead of
%                   Dim by nSubj.
%--------------------------------------------------------------------------
% OUTPUT
% data  an array of size Dim by nSubj as default.
%       (Note that the shape of the array can be changed using
%        the shape of array parameter.)
%--------------------------------------------------------------------------
% EXAMPLES
%--------------------------------------------------------------------------
% AUTHORS:
% Fabian Telschow

% get size of the base_fields and the number of base_fields.
sBase = size(base_fields);
nBase = sBase(end);

% index structure to deal with different dimensions
indexD  = repmat( {':'}, 1, length(sBase)-1 ); 

% Get weights for the multiplier bootstrap
multiplier = normrnd( 0, 1, [ nBase, nSubj ] );

% reshape and and standardize the field, such that it has unit variance
base_fields = reshape( base_fields, prod( sBase(1:end-1) ), nBase );
% normalize the residuals
base_fields = ( base_fields - mean( base_fields, 2 ) ) ...
                                       ./ sqrt( sum( base_fields.^2, 2 ) );

data = zeros( [ sBase( 1:end-1 ), nSubj ] );

for i = 1:nSubj
    % get the bootstrapped process
    data( indexD{:}, i ) = reshape( base_fields * multiplier( :, i ),...
                                        sBase( 1:end-1 ) );
end