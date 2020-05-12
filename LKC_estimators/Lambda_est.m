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
% % 1D 
% nvox = 50; D = 1; FWHM = 3; nsubj = 100;
% lat_data = normrnd(0,1,nvox, nsubj);
% Lambda_est( lat_data, FWHM, D )
%
% % % Stationary evaluation
% FWHM = 10; D = 1; sigma = FWHM2sigma(FWHM); nvox = 100; nsubj = 10000;
% Lambda_theory = diag(repmat(sigma^(-2),1,D))/2; lat_data = normrnd(0,1,nvox,nsubj);
% Lambda_estimates = Lambda_est( lat_data, FWHM, D );
% plot(Lambda_estimates); hold on; abline('h', Lambda_theory);
%
% % 2D
% Dim = [10,10]; D = 2; FWHM = 3; nsubj = 20;
% lat_data = normrnd(0,1,[Dim, nsubj]);
% Lambda_est( lat_data, FWHM, D)
%
% FWHM = 10; sigma = FWHM2sigma(FWHM); Dim = [100,100]; D = length(Dim); nsubj = 500;
% Lambda_theory = diag(repmat(sigma^(-2),1,D))/2; lat_data = normrnd(0,1,[Dim,nsubj]);
% Lambda_estimates = Lambda_est( lat_data, FWHM, D );
% first_entry = Lambda_estimates(1,1,:,:);
% plot(first_entry(:)); hold on; abline('h', Lambda_theory(1,1));
% onetwo_entry = Lambda_estimates(1,2,:,:);
% plot(onetwo_entry(:)) % should be around zero!
%
% % 3D
% FWHM = 6; sigma = FWHM2sigma(FWHM); Dim = [5,5,5]; D = length(Dim); nsubj = 100;
% Lambda_theory = diag(repmat(sigma^(-2),1,D))/2; lat_data = normrnd(0,1,[Dim,nsubj]);
% tic; Lambda_estimates = Lambda_est( lat_data, FWHM, D ); toc
% first_entry = Lambda_estimates(1,1,:,:);
% plot(first_entry(:)); hold on; abline('h', Lambda_theory(1,1));
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if nargin < 4
   resAdd = 1; 
end
if nargin < 5
    h = 0.00001;
end
s_lat_data = size(lat_data);

Dim = s_lat_data(1:end-1);
if length(Dim) ~= D
    error('Incorrect dimension')
end

if resAdd > 1 || resAdd == 0
    spacing = 1/(1+resAdd);
else
    spacing = resAdd;
end

fieldseval = convfield( lat_data, FWHM, spacing, D );
fields_std = sqrt(var(fieldseval, 0, D+1));

fieldseval = fieldseval./fields_std;
clear fields_std

index  = repmat( {':'}, 1, D );     % get variable domain counter

if D == 1
    fieldsplush = convfield( lat_data, FWHM, spacing, D, 0, 1, -1, h );
    fieldsplush_std = sqrt(var(fieldsplush, 0, D+1));
    fieldsplush = fieldsplush./fieldsplush_std;
    clear fieldsplush_std
    
    derivs = (fieldsplush - fieldseval)/h;
    clear fieldsplush fieldseval
    
    Lambda_array = vectcov(derivs, derivs);
    clear derivs
else
    derivs = zeros( [D, s_lat_data] );
    for d = 1:D
        fieldsplush = convfield( lat_data, FWHM, spacing, D, 0, 0, -1, h*sbasis(d,D) );
        fieldsplush_std = sqrt(var(fieldsplush, 0, D+1));
        fieldsplush = fieldsplush./fieldsplush_std;
        clear fieldsplush_std
        derivs( d, index{:}, :) = (fieldsplush - fieldseval)/h;
    end
    clear fieldsplush fieldseval

    Lambda_array = zeros([D,D,Dim]);
    for d1 = 1:D
        for d2 = 1:D
            Lambda_array(d1,d2,index{:}) = vectcov(derivs(d1, index{:}, :), derivs(d2, index{:}, :), D+2);
        end
    end
end

end

