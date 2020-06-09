function divided_mask = mask_highres_divide( mask, divide )
% MASK_HIGHRES_CFIELD( mask, resadd, enlarge ) divides each of the voxels 
% of the original mask into cubes of dimension D and with side length divide
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   mask 
%   resadd
%--------------------------------------------------------------------------
% OUTPUT
%   divided_mask     a high resolution mask where each of the voxels of the
%       original mask have been divided into cubes of dimension D and with
%       side length divide
%--------------------------------------------------------------------------
% EXAMPLES
% mask = ones(3); mask(2,1) = 0; divide = 2;
% mask_hr = mask_highres_divide(mask, divide);
% imagesc(mask_hr)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Get dimensions
Dim = size(mask); D = length(Dim);

%%  add/check optional values
%--------------------------------------------------------------------------
% if ~exist( 'enlarge', 'var' )
%    % default option of enlarge
%    enlarge = ceil(resadd/2);
% end

%%  Main Function Loop
%--------------------------------------------------------------------------
dividedDim = divide*Dim;
index = cell( [ 1 D ] );
for d = 1:D
    index{d} = divide:divide:dividedDim(d);
end
D_vec = divide*ones(1,D);
divided_mask = zeros(dividedDim);
divided_mask(index{:}) = mask;
divided_mask = convn(divided_mask, ones(D_vec));

inner_index = cell( [ 1 D ] );
for d = 1:D
    inner_index{d} = divide:(dividedDim(d)+divide-1);
end
divided_mask = divided_mask(inner_index{:});

% D_vec = D*ones(1,D);
% mask_hr = convn(mask, ones(D_vec));

% Get the box kernel
% RK = @(y)boxker(y, 0.5, 1);
% spacing = 1/(1+resadd);
% 
% % Need to change this once have added a non_separable option to convfield
% % Convolve the mask with a the box kernel to return the high resolution
% % version
% mask_hr = convfield_dep(mask, RK, spacing, D, 0, 0, 1);

end

