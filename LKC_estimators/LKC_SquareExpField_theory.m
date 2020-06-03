function [ LKC, resel ] = LKC_SquareExpField_theory( FWHM, dim, type )
% This function computes the theoretical LKCs for the Gaussian random field
% with the square exponential covariance function
% C(s) = exp( -s^2 / 2 / FWHM^2 )
% over an rectangle S = T_1 x ... x T_D for D < 4.
% In other words it computes the LKCs of smoothed white noise.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   FWHM  a numeric giving the FWHM for the Gaussian smoothing kernel
%         used to smooth white noise.
%   dim   a 1 by D vector giving the length of the sides of the
%         rectangle in each dimension. Note that length(dim) < 4 is
%         required.
%
% Optional
%   type   a string indicating whether FWHM or sigma should be usd.
%          Options are 'sigma', 'fwhm'.
%          Default 'fwhm'.
%--------------------------------------------------------------------------
% OUTPUT
%   LKC    a 1 by D vector containing the true LKCs computed from the
%          analytical formula.
%   resel  a 1 by D vector of the LKCs transformed into Worsley's resels
%--------------------------------------------------------------------------
% DEVELOPER TODOs:
% -------------------------------------------------------------------------
% EXAMPLES
%--------------------------------------------------------------------------
% AUTHORS: Fabian Telschow, Wenyi Lin
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Check mandatory input and get constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check whether the input is valid
if length( dim ) > 3
    error( "The function supports only rectangles up to dimension 3." )
end

% compute the dimension
D = length( dim ); 
if D == 2 && ( dim(1) == 1 || dim(2) == 1 )
    D = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% add/check optional values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist( 'type', 'var' )
   % default number of resolution increasing voxels between observed voxels
   type = "fwhm";
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% main function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theoretical value for Lambda matrix. Note that because of
% stationarity it is constant times identity accross the domain.
switch type
    case "sigma"
       alpha =  1 / ( 4 * FWHM^2 );
    case 'fwhm'
       alpha = 1 / ( 4 * FWHM2sigma( FWHM )^2 );
end

% compute the LKCs by multiplication with certain volumes compare:
% Worsley et al (2004) Unified univariate and multivariate random field
% theory
switch D
   case 1
       LKC = dim(1) * sqrt( 2 * alpha );
   case 2
       LKC = [ sum( dim ); prod( dim ) ] .* ...
                            [ sqrt( 2 * alpha ); 2 * alpha ];
   case 3
       LKC = [ sum( dim ) / 2; dim(1)*dim(2) + dim(1)*dim(3)...
                    + dim(2)*dim(3); prod( dim ) ] .* ...
                      [ sqrt( 2 * alpha ); 2 * alpha; ( 2 * alpha )^(3/2) ];
end
LKC = LKC';

% compute the resels
resel = ( LKC ./ sqrt( ( 4 * log(2) ).^(1:D) ) );

return