function LD = LDcalc( lat_data, FWHM, D, mask, spacing, h )
% LDcalc( lat_data, FWHM, mask ) calculates the Dth LKC for a D-dimensional
% convolution t field generated from lattice data with an isotropic Gaussian
% kernel with given FWHM.
%--------------------------------------------------------------------------
% ARGUMENTS
% lat_data      a Dim by nsubj array corresponding to the lattice data
% FWHM          the FWHM of the smoothing kernel
% mask          a 0/1 Dim by nsubj array with 1s indicating the mask
%--------------------------------------------------------------------------
% OUTPUT
% LD            an estimate of the Dth LKC
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if nargin < 5
   spacing = 1; 
end
if nargin < 6
    h = 0.00001;
end
s_lat_data = size(lat_data);

Dim = s_lat_data(1:end-1);
if length(Dim) ~= D
    error('Incorrect dimension')
end
if nargin < 4
    if D == 1
        mask = ones([Dim,1]);
    else
        mask = ones(Dim);
    end
elseif (D == 1) && nargin == 4
    if size(mask, 1) == 1
        mask = mask';
    end
end

fieldseval = convfield( lat_data, FWHM, spacing, D );
fields_std = sqrt(var(fieldseval, 0, D+1));

fieldseval = fieldseval./fields_std;
clear fields_std

index  = repmat( {':'}, 1, D+1 );     % get variable domain counter

if D == 1
    fieldsplush = convfield( lat_data, FWHM, spacing, D, 0, 1, h );
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
        fieldsplush = convfield( lat_data, FWHM, spacing, D, 0, 1, h*sbasis(d,D) );
        fieldsplush_std = sqrt(var(fieldsplush, 0, D+1));
        fieldsplush = fieldsplush./fieldsplush_std;
        clear fieldsplush_std
        derivs( d, index{:}) = (fieldsplush - fieldseval)/h;
    end
    clear fieldsplush fieldseval

    Lambda_array = zeros(D,D,prod(Dim));
    for d1 = 1:D
        for d2 = 1:D
            Lambda_array(d1,d2) = vectcov(derivs(d1, index{:}), derivs(d1, index{:}), 5);
        end
    end
end

if D == 1
    Lambda_array = Lambda_array';
    mask = mask';
end

LD = sum((abs(vectdet(Lambda_array)).^(1/2)))*spacing;  % Need to deal with the masking 

% LD = sum((abs(vectdet(Lambda_array)).^(1/2)).*mask)*spacing; 

% No scaling needed above as each voxel has volume 1

end

% if D == 1
%     Dim = s_lat_data(1:end-1);
%     if Dim(1) ~= 1 && Dim(2) ~= 1 || length(Dim) > 2
%         error('Incorrect dimension')
%     end
% else
%     Dim = s_lat_data(2:end-1);
%     if length(Dim) ~= D
%         error('Incorrect dimension')
%     end
% end