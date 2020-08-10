function data = multiplier_field( base_fields, nsubj )
% multiplier_field( base_fields, nsubj ) generates samples from a mean zero
% unit variance field given by the base_fields using Gaussian multipliers.
%--------------------------------------------------------------------------
% ARGUMENTS
% base_fields   an object of class Field containing base_fields.fibersize
%               base_fields.
% nsubj         The number of subjects to be drawn from the base_fields
%               gaussian multiplier.
%--------------------------------------------------------------------------
% OUTPUT
% data   an object of class Field containing the nsubj sample fields
%--------------------------------------------------------------------------
% EXAMPLES
%--------------------------------------------------------------------------
% AUTHORS:
% Fabian Telschow


data2 = base_fields;
base_fields = data2.field;
% get size of the base_fields and the number of base_fields.
sBase = size(base_fields);
nBase = sBase(end);

% index structure to deal with different dimensions
indexD  = repmat( {':'}, 1, length(sBase)-1 );
D = length(sBase)-1;

% Get weights for the multiplier bootstrap
multiplier = normrnd( 0, 1, [ nBase, nsubj ] );

% reshape and and standardize the field, such that it has unit variance
base_fields = reshape( base_fields, prod( sBase(1:end-1) ), nBase );
% normalize the residuals
base_fields = ( base_fields - mean( base_fields, 2 ) ) ...
                                       ./ sqrt( sum( base_fields.^2, 2 ) );

data = zeros( [ sBase( 1:end-1 ), nsubj ] );

for i = 1:nsubj
    % get the bootstrapped process
    if D>1
        data( indexD{:}, i ) = reshape( base_fields * multiplier( :, i ),...
                                        sBase( 1:end-1 ) );
    else
        data( indexD{:}, i ) = base_fields * multiplier( :, i );
    end
end

data2.field = data;
data = data2;