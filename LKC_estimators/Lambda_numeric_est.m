function [Lambda_array, xvals ] = Lambda_numeric_est( lat_data, Kernel,...
                                                      resadd, enlarge, h )
% LAMBDA_NUMERIC_EST( lat_data, FWHM, resadd, enlarge, h ) calculates an 
% estimate of Lambda(v) = cov(\nabla X(v)) at each voxel
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   lat_data data array T_1 x ... x T_D x N. Last index enumerates the
%            number of subjects. Note that N > 1 is required!
%   Kernel   array 1x1 or 1xD containing the FWHM for different directions
%            for smoothing with a Gaussian kernel, if numeric an isotropic
%            kernel is assumed.
% Optional
%   resadd   integer denoting the amount of voxels padded between existing
%            voxels to increase resolution
%   enlarge  integer denoting the amount of voxels padded between existing
%            voxels to increase resolution
%   h        the h used to calculate the derivatives i.e. via
%            (X(v+h)-X(v))/h. Default is 0.00001. Avoid taking h to be
%            too small for numerical precision reasons
%--------------------------------------------------------------------------
% OUTPUT
% Lambda_array  a D by D by prod(Dim) array giving the Lambda matrix at each
%               voxel
% -------------------------------------------------------------------------
% DEVELOPER TODOs:
%--------------------------------------------------------------------------
% EXAMPLES
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport, Fabian Telschow
%--------------------------------------------------------------------------

%% Check input and get important constants from the mandatory input
%--------------------------------------------------------------------------
% Size of the domain
s_lat_data = size( lat_data );

% Dimension of the domain
D = length( s_lat_data( 1:end-1 ) );

% Get variable domain counter
index  = repmat( {':'}, 1, D );

% Check that method is implemented for dimension D
if D > 3
    error( 'D must be < 4. Higher dimensional domains have not been implemented')
end

%% add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'resadd', 'var' )
   % default number of resolution increasing voxels between observed voxels
   resadd = 1;
end

if ~exist( 'enlarge', 'var' )
   % default number of resolution increasing voxels between observed voxels
   enlarge = ceil( resadd / 2 );
end

if  ~exist( 'h', 'var' )
    h = 0.0000001; % Note don't take h to be substantially lower than this 
                 % or there will be numerical precision issues
end


%% main function
%--------------------------------------------------------------------------
% Evaluate the value of the convolution field
[ fieldseval, xvals, Kernel ] = convfield( lat_data, Kernel, resadd,...
                                                           D, 0, enlarge ); 
% Get the standard deviation of the fields everywhere                                 
fields_std = sqrt( var( fieldseval, 0, D+1 ) );

% Obtain a variance 1 field
fieldseval = fieldseval ./ fields_std; 
clear fields_std

% Get the size of the high resolution field
Dimhr = size( fieldseval );

% Main loop
if D == 1
    % Define the shifted kernel
    Kernelh = Kernel;
    Kernelh.adjust = h;

    % Obtain the variance 1 field at an offset of h everywhere
    fieldsplush = convfield( lat_data, Kernelh, resadd, D, 0, enlarge );
    fieldsplush_std = sqrt( var( fieldsplush, 0, D+1 ) );
    fieldsplush     = fieldsplush ./ fieldsplush_std;
    clear fieldsplush_std
    
    % Calculate the derivatives of the variance 1 fields
    derivs = ( fieldsplush - fieldseval ) / h;
    clear fieldsplush fieldseval
    
    % Calculate Lambda
    Lambda_array = vectcov( derivs, derivs, 2, 0 );
    clear derivs
    
elseif D == 2 || D == 3
    % allocate variable for the derivatives of the fields
    derivs = zeros( [ Dimhr, D ] );
    
    for d = 1:D
        % define the Kernel object to be for shifted h along an offset of
        % the d-th standard direction or R^D
        Kernel.adjust = h*sbasis( d, D );

        % Obtain the variance 1 field at an offset of h*e_d everywhere
        fieldsplush = convfield( lat_data, Kernel, resadd, D, 0, enlarge );                                          
        fieldsplush_std = sqrt( var( fieldsplush, 0, D+1 ) );
        fieldsplush = fieldsplush ./ fieldsplush_std;
        clear fieldsplush_std
        
        % Calculate the partial derivatives in the e_d th direction
        derivs( index{:}, :, d ) = ( fieldsplush - fieldseval ) / h;
    end
    clear fieldsplush fieldseval

    % Calculate Lambda
    Lambda_array = zeros( [ Dimhr(1:end-1),  D, D ] );
    for d1 = 1:D
        for d2 = 1:D
            % Calculate the covariance between the d1th partial derivative
            % and the d2th partial derivative
            Lambda_array( index{:}, d1, d2 ) = vectcov( ...
                   derivs(index{:}, :, d1), derivs(index{:}, :, d2 ),D+1);
        end
    end
else
    error('Not working in D > 3 yet')
end

end