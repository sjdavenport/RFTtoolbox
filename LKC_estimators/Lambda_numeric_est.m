function Lambda_array = Lambda_numeric_est( lat_data, Kernel, resAdd,...
                                                  enlarge, h )
% LAMBDA_EST( lat_data, FWHM, D, spacing, h ) calculates an estimate of 
% Lambda(v) = cov(\nabla X(v)) at each voxel
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   lat_data data array T_1 x ... x T_D x N. Last index enumerates the
%            samples. Note that N > 1 is required!
%   Kernel   array 1x1 or 1xD containing the FWHM for different directions
%            for smoothing with a Gaussian kernel, if numeric an isotropic
%            kernel is assumed.
% Optional
%   resAdd   integer denoting the amount of voxels padded between existing
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
% % 1D (stationary example)
% FWHM = 10; D = 1; sigma = FWHM2sigma(FWHM); nvox = 100; nsubj = 10000;
% Lambda_theory = diag(repmat(sigma^(-2),1,D))/2; lat_data = normrnd(0,1,nvox,nsubj);
% Lambda_estimates = Lambda_est( lat_data, FWHM, D );
% Lambda_theory
% mean(Lambda_estimates(5:95)) % Not the ends because of edge effects
% plot(Lambda_estimates); hold on; abline('h', Lambda_theory);
%
% % 2D (stationary example)
% FWHM = 5; sigma = FWHM2sigma(FWHM); Dim = [100,100]; D = length(Dim); nsubj = 50;
% Lambda_theory = diag(repmat(sigma^(-2),1,D))/2; lat_data = normrnd(0,1,[Dim,nsubj]);
% Lambda_estimates = Lambda_est( lat_data, FWHM, D );
% first_entry = Lambda_estimates(1,1,5:95,5:95);
% Lambda_theory(1,1)
% mean(first_entry(:))
% onetwo_entry = Lambda_estimates(1,2,5:95,5:95);
% Lambda_theory(1,2)
% mean(onetwo_entry(:))
% plot(first_entry(:)); hold on; abline('h', Lambda_theory(1,1));
%
% % 3D (stationary example)
% FWHM = 10; sigma = FWHM2sigma(FWHM); Dim = [50,50,50]; D = length(Dim); nsubj = 50;
% Lambda_theory = diag(repmat(sigma^(-2),1,D))/2; lat_data = normrnd(0,1,[Dim,nsubj]);
% Lambda_estimates = Lambda_est( lat_data, FWHM, D );
% first_entry = Lambda_estimates(1,1,5:45,5:45,5:45);
% Lambda_theory(1,1)
% mean(first_entry(:))
% onetwo_entry = Lambda_estimates(1,2,5:45,5:45,5:45);
% Lambda_theory(1,2)
% mean(onetwo_entry(:))
% plot(first_entry(:)); hold on; abline('h', Lambda_theory(1,1));
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport, Fabian Telschow
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Check input and get important constants from the mandatory input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get constants from the mandatory input
% size of the domain
s_lat_data = size( lat_data );

% dimension of the domain
D = length( s_lat_data( 1:end-1 ) );

% size of domain
Dim = s_lat_data( 1:end-1 );

% get number of subjects/samples
nsubj = s_lat_data( D + 1 );

% get variable domain counter
index  = repmat( {':'}, 1, D );

% check that method is implemented for dimension D
if D > 3
    error( 'D must be < 4. Higher dimensional domains have not been implemented')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% add/check optional values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist( 'resAdd', 'var' )
   % default number of resolution increasing voxels between observed voxels
   resAdd = 1;
end

if ~exist( 'enlarge', 'var' )
   % default number of resolution increasing voxels between observed voxels
   enlarge = ceil( resAdd / 2 );
end

if  ~exist( 'h', 'var' )
    h = 0.00001; % Note don't take h to be substantially lower than this 
                 % or there will be numerical precision issues
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% main function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate the value of the convolution field
[ fieldseval, ~, Kernel ] = convfield_struct( lat_data, Kernel, resAdd,...
                                              D, 0, enlarge ); 
% Get the standard derivation of the fields everywhere                                 
fields_std = sqrt( var( fieldseval, 0, D+1 ) );

% Obtain a variance 1 field
fieldseval = fieldseval ./ fields_std; 
clear fields_std

% Define the shifted kernel
Kernelh = struct();
Kernelh.adjust_kernel = ;



% Main loop
if D == 1
    % Obtain the variance 1 field at an offset of h everywhere
    fieldsplush = convfield_struct( lat_data, Kernelh, resAdd, D, 0,...
                                     enlarge );
    fieldsplush_std = sqrt(var(fieldsplush, 0, D+1));
    fieldsplush = fieldsplush./fieldsplush_std;
    clear fieldsplush_std
    
    % Calculate the derivatives of the variance 1 fields
    derivs = (fieldsplush - fieldseval)/h;
    clear fieldsplush fieldseval
    
    % Calculate Lambda
    Lambda_array = vectcov(derivs, derivs);
    clear derivs
elseif D == 2 || D == 3
    derivs = zeros( [D, s_lat_data] );
    for d = 1:D
        % Obtain the variance 1 field at an offset of h*e_d everywhere
        % where e_d is the standard basis vector in the dth direction
        fieldsplush = convfield( lat_data, FWHM, spacing, D, 0, -1, h*sbasis(d,D) );
        fieldsplush_std = sqrt(var(fieldsplush, 0, D+1));
        fieldsplush = fieldsplush./fieldsplush_std;
        clear fieldsplush_std
        
        % Calculate the partial derivatives in the e_d th direction
        derivs( d, index{:}, :) = (fieldsplush - fieldseval)/h;
    end
    clear fieldsplush fieldseval

    % Calculate Lambda
    Lambda_array = zeros([D,D,Dim]);
    for d1 = 1:D
        for d2 = 1:D
            % Calculate the covariance between the d1th partial derivative
            % and the d2th partial derivative
            Lambda_array(d1,d2,index{:}) = vectcov(derivs(d1, index{:}, :), derivs(d2, index{:}, :), D+2);
        end
    end
else
    error('Not working in D > 3 yet')
end

end

