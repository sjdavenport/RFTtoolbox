function Lambda_array = Lambda_est( lat_data, FWHM, D, resAdd, h )
% LAMBDA_EST( lat_data, FWHM, D, spacing, h ) calculates an estimate of 
% Lambda(v) = cov(\nabla X(v)) at each voxel
%--------------------------------------------------------------------------
% ARGUMENTS
% lat_data      a Dim by nsubj array corresponding to the lattice data.
%               note must have nsubj > 1
% FWHM          the FWHM of the smoothing kernel
% D             the dimension
% resAdd        the resolution increase between the voxels. I.e 0 adds no
%               point in between, 1 adds one etc
% h             the h used to calculate the derivatives i.e. via
%               (X(v+h)-X(v))/h. Default is 0.00001. Avoid taking h to be
%               too small for numerical precision reasons
%--------------------------------------------------------------------------
% OUTPUT
% Lambda_array  a D by D by prod(Dim) array giving the Lambda matrix at each
%               voxel
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
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if nargin < 4
   resAdd = 1; %Default is a resolution of 1, i.e. a spacing of 0.5
end
if nargin < 5
    h = 0.00001; % Note don't take h to be substantially lower than this 
                 %or there will be numerical precision issues
end
s_lat_data = size(lat_data);

Dim = s_lat_data(1:end-1);
if length(Dim) ~= D
    error('Incorrect dimension')
end

% Define spacing in terms of the resolution
if resAdd > 1 || resAdd == 0
    spacing = 1/(1+resAdd);
else
    spacing = resAdd;
end

% Evaluate the value of the convolution field
fieldseval1 = convfield( lat_data, FWHM, spacing, D );
fields_std = sqrt(var(fieldseval1, 0, D+1)); % Get the standard derivation of the fields everywhere

%Obtain a variance 1 field
fieldseval1 = fieldseval1./fields_std; 
clear fields_std

% Get variable domain counter
index  = repmat( {':'}, 1, D );     

% Main loop
if D == 1
    % Obtain the variance 1 field at an offset of h everywhere
    fieldsplush = convfield( lat_data, FWHM, spacing, D, 0, -1, h );
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

