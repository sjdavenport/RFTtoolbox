function mfield = multiplier_field( base_fields, nsubj, normalize, multiplier )
% multiplier_field( base_fields, nsubj ) generates samples from a mean zero
% unit variance field given by the base_fields using Gaussian multipliers.
%--------------------------------------------------------------------------
% ARGUMENTS
% base_fields   an object of class Field containing base_fields.fibersize
%               base_fields.
% nsubj         The number of subjects to be drawn from the base_fields
%               gaussian multiplier.
% normalize     a boolean, if TRUE the field is normalized to
%               have variance one. Default TRUE.
% multiplier    
%--------------------------------------------------------------------------
% OUTPUT
% data   an object of class Field containing the nsubj sample fields
%--------------------------------------------------------------------------
% EXAMPLES
% % Grid for the domain
% x = (0:50) / 50;
% y = (0:50) / 50;
% [X, Y] = meshgrid(x, y);
% 
% % Mask
% mask = true(length(x), length(y));
% 
% % Subjects
% nsubj = 1e3;
% 
% % Array containing polynomials up to order 2
% Polynomials = cat(3, ones(length(x), length(y)), X, Y, X.*Y, X.*X, Y.*Y);
% 
% % Generate a Basis functions as polynomials 
% PolyBase = Field(Polynomials, mask);
% 
% % Generate data with variance 1
% mfields_const = multiplier_field(PolyBase, nsubj, true, "gaussian");
% imagesc(mfields_const(:, :, 100));
% 
% % Plot the variance of the field
% imagesc(std(mfields_const)); colorbar;
% 
% % Generate data with non constant variance given by the basis functions
% mfields_nonconst = multiplier_field(PolyBase, 20, false, "gaussian");
% imagesc(mfields_nonconst(:, :, 20))
% 
% % Plot the variance of the field
% var(mfields_nonconst)
% 
% % Image basis
% base_fields = Field(randn([91,109,40]), true([91,109]));
% multiplier_field( base_fields, 10 ) 
%--------------------------------------------------------------------------
% AUTHORS:
% Fabian Telschow
%--------------------------------------------------------------------------

%% Check mandatory input
%--------------------------------------------------------------------------

%% Check optional input
%--------------------------------------------------------------------------
% Fill the default for the normaliez option
if ~exist( 'normalize', 'var')
    normalize = true;
end

% Fill the default for the normalize option
if ~exist( 'multiplier', 'var')
    multiplier = "gaussian";
end

%% Main function
%--------------------------------------------------------------------------
% save the field paramters from the basis field into a temporary variable
base_fiber = base_fields.field;

% get size of the base_fields and the number of base_fields.
sBase = base_fields.masksize;
nBase = base_fields.fibersize;

% Get weights for the multiplier bootstrap (atm only Gaussian implemented)
if( multiplier == "gaussian" )
    multiplier = normrnd( 0, 1, [ nBase, nsubj ] );
end

% reshape and and standardize the field, such that it has unit variance
base_fiber = reshape( base_fiber, prod( sBase ), nBase );

% normalize the basefields to yield variance one
if( normalize )
    base_fiber = base_fiber ./ sqrt( sum( base_fiber.^2, 2 ) );
end

% Get the random multiplier fields
mfield = base_fields;
mfield.field = reshape( base_fiber * multiplier, [ sBase, nsubj ] );

return