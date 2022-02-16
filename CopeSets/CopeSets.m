function [ thresh, quantiles, hatmu, hatsigma, len_bdry ] = CopeSets( field, c, lvls, quantEstim,...
                                                        bdry_type, center, normalize, delta )
% CopeSets( field, c, lvls, quantEstim, bdry_type, center, normalize,
% delta ) computes CoPe sets for the mean of a random field.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  field  random field over a domain in R^D, it is an (D+1)-dimensional
%         array, where the last dimension enumerates the samples
%  c      threshold for excursions
%  lvls   vector containing the required confidence levels. Must be
%         between 0 and 1.
% Optional
%  quantEstim  structure containing the name and the parameters for the
%              quantile estimation method. Choices:
%               {
%                quantEstim.name = 'multiplierbootstrap'
%                quantEstim.params:
%                   Mboot:     amount of bootstrap replicates (default=5e3)
%                   method:    option for the bootstrap estimator (default='t')
%               }
%  bdry_type  currently 'linear', 'eroddilate' or 'true' are supported.
%             Moreover, we have the 'interval'
%  center     option to center the field using the sample mean (default=1)
%  normalize  option to normalize the field by sample variance (default=1)
%  delta      required, if bdry_type is equal to 'true'. This is the true
%             population mean function given on a D-dimensional array
%--------------------------------------------------------------------------
% OUTPUT
%  thresh    the threshold lower and upper for the sample mean in order to
%            be in the estimated lower and upper excursion sets 
%  quantile  the bootstrapped quantile of the maximum distribution of the 
%            input processes
%  hatmu  the sample mean of the fields
%  hatsigma  the sample variance of the fields
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR:  Fabian Telschow
%--------------------------------------------------------------------------

%% Check mandatory input
%--------------------------------------------------------------------------
% Check confidence levels input
if any( lvls >= 1 ) || any( lvls <= 0 )
    error( "The vector lvls needs to have entries between 0 and 1!" )
end

% Fill in unset optional values.
switch nargin
    case 3
        quantEstim = struct( 'name', "MultiplierBootstrap",...
                             'params', struct( 'Mboot', 5e3,...
                                               'weights', "rademacher",...   
                                               'method', 't' )...
                                              );
        bdry_type = 'linear';
        center    = 1;
        normalize = 1;
        delta     = 666;
    case 4
        bdry_type = 'linear';
        center    = 1;
        normalize = 1;
        delta     = 666;
    case 5
        center    = 1;
        normalize = 1;
        delta     = 666;
    case 6
        normalize = 1;
        delta     = 666;
    case 7
        delta     = 666;
end

% Get useful constants
D     = field.D;
N     = field.fibersize;
index = repmat( {':'}, 1, length( field.fieldsize ) );
sF    = field.masksize;


%% Main function
%--------------------------------------------------------------------------
% Compute mean curve and variance
hatmu = mean( field );
hatsigma = std( field );

%%% Compute the process on the boundary and its mask
switch bdry_type
    case 'linear'
        F_bdry = BdryEst_linear( field, c );
        mask   = ones( [ size( F_bdry, 1 ) 1 ] );
    case 'erodilation'
        F_bdry = BdryEst_erodedilate( field, c );
        mask   = ones( [ size( F_bdry, 1 ) 1 ] );
    case 'true'
        F_bdry = BdryEst_linear( field, c, delta );
        mask   = ones( [ size( F_bdry, 1 ) 1 ] );
end

% Compute length of estimated boundary
len_bdry = sum( mask );

%%% convert F_bdry to residuals if neccessary
% Center/normalize, if required
if center
    F_bdry = F_bdry - mean( F_bdry, 2 );
end
if normalize
   F_bdry = F_bdry ./ std( F_bdry, 0, 2 );
end

%%% Estimate the quantiles of the Gaussian process on the boundary
if strcmp( quantEstim.name, 'MultiplierBootstrap' )
    quantiles = MultiplierBoots( F_bdry,...
                                 lvls,...
                                 quantEstim.params.Mboot,...
                                 mask,...
                                 quantEstim.params.weights,...
                                 quantEstim.params.method );
else
    error( "Please specify a valid method for quantile estimation." )
end

%%% Compute upper and lower threshold for CoPE sets
thresh = ones( [ sF( 1:end ) length( lvls ) 2 ] );
% Lower threshold
thresh( index{:}, 1 ) = c - repmat( shiftdim(  quantiles,...
                                              -D + 1 ),...
                                    [ sF 1 ] )...
                        .* repmat( hatsigma.field,...
                                   [ 1 1 length( lvls ) ] ) / sqrt( N );
% Upper threshold
thresh( index{:}, 2 ) = c + repmat( shiftdim(  quantiles,...
                                              -D+1),...
                                              [ sF 1 ] )...
                        .* repmat( hatsigma.field,...
                                   [ 1 1 length( lvls ) ] ) / sqrt( N );