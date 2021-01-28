function [ quantiles, bootMax ] = MultiplierBoots( R, alpha, Mboot, mask, weights, method )
% MultiplierBoots( R, alpha, Mboot, mask, weights, method ) estimates the
% alpha-quantiles of the maximum of a random field using the multiplier
% bootstrap.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   R        random field over a domain in R^D, it is an (D+1)-dimensional
%            array, where the last dimension enumerates the samples
%   alpha    vector of quantiles
% Optional
%   Mboot    amount of bootstrap replicates. Default = 5e3.
%   mask     option of specifying an boolean array of size of the first D
%            components of size( R ) containing the locations used for the
%            multplier bootstrap.
%            Default = ones( size(R)[ 1:D ] )
%  weights   option type of weights used as multipliers. Options are
%               - 'gaussian': Gaussian random variables 
%               - 'rademacher': Rademacher random variables ( Default )
%  method:   options are 't' for normalizing by the bootstrapped variance
%            or 'regular' for not normalizing.
%            Default = 't'.
%
%--------------------------------------------------------------------------
% OUTPUT
%  quantile  the bootstrapped quantile of the maximum distribution of the 
%            input processes
%  bootMax   the bootstrap distribution of the maximum of the process
%
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR:  Fabian Telschow
%--------------------------------------------------------------------------
%% Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Get parameters for simulation
dimR = size( R );
N    = dimR( end );


%% Add/check optional values
%--------------------------------------------------------------------------
% Fill in unset optional values.
switch nargin
    case 2
        Mboot     = 5e3;
        if( length( dimR ) == 2 )
           mask   = ones( [ dimR( 1 ) 1 ] );
        else
           mask   = ones( dimR( 1:end - 1 ) );
        end
        weights   = 'rademacher';
        method    = 't';
        
    case 3
        if( length( dimR ) == 2 )
           mask   = ones( [ dimR(1) 1 ] );
        else
           mask   = ones( dimR( 1:end-1 ) );
        end
        weights   = 'rademacher';
        method    = 't';
        
    case 4
        weights   = 'rademacher';
        method    = 't';
        
    case 5
        method    = 't';
        
end


%%  Main Function
%--------------------------------------------------------------------------
% Combine values of R in the mask to an matrix for faster computation
R       = reshape( R, prod( dimR( 1:end-1 ) ), N );
mask    = repmat( reshape( logical( mask ), prod( dimR( 1:end-1 ) ), 1 ),...
                  [ 1 N ] );

R = R( mask );
R = reshape( R, [ length( R ) / N N ] );        

%%% compute bootstrap replicates
% compute multipliers
if strcmp( weights, 'gaussian' )
    multiplier = normrnd( 0, 1, [ N, Mboot ] );
elseif strcmp( weights, 'rademacher' )
    multiplier = randi( 2, [ N, Mboot ] ) * 2 - 3;
else
    error( 'Error: Please, choose weights from the available options "gaussian" or "rademacher"!' )
end

% compute bootstrapped means
bootMeans      = R * multiplier / N;

if method == 't'
    bootSecMoments = R.^2 * multiplier.^2 / N;
    % we put an abs here to make sure that no NaNs are produced due to machine precision error.
    bootSigma = sqrt( abs( bootSecMoments - bootMeans.^2 ) / ( N - 1 ) * N );
else
    bootSigma = 1;
end

% compute bootstrapped values of maximum
bootMax = max( abs( sqrt( N ) * bootMeans ./ bootSigma ) );

% compute quantiles from the bootstrapp distribution of maximum
quantiles = quantile( bootMax, alpha / 2 );

% bootMax          = max(sqrt(N)*bootMeans./bootSigma);
% quantilesAsym    = [-666 -666];
% quantilesAsym(1) = quantile( bootMax, alpha );
% bootMax          = min(sqrt(N)*bootMeans./bootSigma);
% quantilesAsym(2) = quantile( bootMax, 1-alpha );

