function obj = GaussKernel( obj, FWHM )
% GAUSSKERNEL( obj, FWHM ) constructs from a basic Kernel class
% object a Kernel class object representing a seperable Gaussian
% Kernel and fills all properties.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  obj   an SepKernel object
%  FWHM  either a numeric or a numeric vector of length obj.D
%
%--------------------------------------------------------------------------
% OUTPUT
% obj  an object of class SepKernel representing a seperable
%      Gaussian Kernel.
%
%--------------------------------------------------------------------------
% EXAMPLES
% % create an object of class SepKernel representing an isotropic
% % Gaussian kernel
% D = 2
% sepK = SepKernel( D )
% sepK = GaussKernel( sepK, 6 )
%
% % create an object of class SepKernel representing an seperable
% % Gaussian kernel
% D = 2
% sepK = SepKernel( D )
% sepK = GaussKernel( sepK, [ 6, 2 ] )
%
%--------------------------------------------------------------------------
% Author: Fabian Telschow
%--------------------------------------------------------------------------

%% Check input
%--------------------------------------------------------------------------
% Check input argument and make it a vector
if length( FWHM ) == 1
  FWHM = FWHM * ones( [ 1 obj.D ] );
end

%% Main function
%--------------------------------------------------------------------------

if  length( FWHM ) == obj.D
% Fill the field kernel with Gaussians
for d = 1:obj.D
    obj.kernel{d} = @(x) Gker( x, FWHM(d) );
end

% Fill the field dkernel with derivatives of the Gaussians
for d = 1:obj.D
    obj.dkernel{d} = @(x) Gkerderiv( x, FWHM(d) );
    obj.d2kernel{d} = @(x) Gkerderiv2( x, FWHM(d) );
end

% Fill the field truncation and dtruncation
obj.truncation   = ceil( 3 * FWHM2sigma( FWHM ) );
obj.dtruncation  = obj.truncation;
obj.d2truncation = obj.truncation;
else
  error( "FWHM needs to be either of length 1 or D." )
end

return