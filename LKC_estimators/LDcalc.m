function LD = LDcalc( lat_data, FWHM, D, mask, spacing, h )
% LDcalc( lat_data, FWHM, mask ) calculates the Dth LKC for a D-dimensional
% convolution t field generated from lattice data  with an isotropic Gaussian
% kernel with given FWHM.
%--------------------------------------------------------------------------
% ARGUMENTS
% lat_data      a Dim by nsubj array corresponding to the lattice data
% FWHM          the FWHM of the smoothing kernel
% D             the dimension
% mask          a 0/1 Dim by nsubj array with 1s indicating the mask
% spacing       the spacing between the voxels
% h             the h used to calculate the deratives i.e. via
%               (X(v+h)-X(v))/h. Default is 0.00001. Avoid taking h to be
%               too small for numerical precision reasons
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

Lambda_array = Lambda_est( lat_data, FWHM, D, resAdd, h );

if D == 1
    Lambda_array = Lambda_array';
    mask = mask';
end

LD = sum((abs(vectdet(Lambda_array)).^(1/2)))*spacing;  % Need to deal with the masking 

% smaller mask version
if D == 1
    
elseif D == 2
    
elseif D == 3
    
end

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