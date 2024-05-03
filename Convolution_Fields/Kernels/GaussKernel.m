function obj = GaussKernel(D, FWHM, adjust)
% GAUSSKERNEL( D, FWHM, adjust ) constructs from a basic Kernel class
% object a Kernel class object representing a seperable Gaussian
% Kernel and fills all properties.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  D     the dimension of the Gauss kernel
%  FWHM  either a numeric or a numeric vector of length D
% Optional
%  adjust fills the adjust field of the SepKernel class
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
% sepK = GaussKernel( D, 6 )
%
% % create an object of class SepKernel representing an seperable
% % Gaussian kernel
% D = 2
% sepK = GaussKernel( D, [ 6, 2 ] )
%
%--------------------------------------------------------------------------
% Author: Fabian Telschow
%--------------------------------------------------------------------------

%% Check input
%--------------------------------------------------------------------------
% Check input argument and make it a vector
if length( FWHM ) == 1
  FWHM = FWHM * ones( [ 1 D ] );
end

if ~exist('adjust', 'var')
    adjust = zeros([1 D]);
end

%% Main function
%--------------------------------------------------------------------------
% Create a SepKernel object to be filled with the Gaussian kernels
obj = SepKernel(D);

if  length( FWHM ) == D
    % Fill the field kernel with Gaussians
    for d = 1:D
        obj.kernel{d} = @(x) Gker( x, FWHM(d) );
    end

    % Fill the field dkernel with derivatives of the Gaussians
    for d = 1:D
        obj.dkernel{d}  = @(x) Gkerderiv( x, FWHM(d) );
        obj.d2kernel{d} = @(x) Gkerderiv2( x, FWHM(d) );
    end

    % Fill the field truncation and dtruncation
    obj.truncation   = ceil( 4 * FWHM2sigma( FWHM ) );
    obj.dtruncation  = obj.truncation;
    obj.d2truncation = obj.truncation;

    % Fill the adjust field
    obj.adjust = adjust;

else
  error( "FWHM needs to be either of length 1 or D." )
end

return