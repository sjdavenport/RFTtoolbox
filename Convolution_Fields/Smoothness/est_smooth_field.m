function [ fwhm_est_forman, fwhm_est_kiebel, Lambda_est, sigma_est] = est_smooth_field( field, df )
% EST_SMOOTH estimates the smoothness of a process. NEED TO DO THE
% CROSS TERMS!! Would be good to study exactly how biased this is.
%--------------------------------------------------------------------------
% ARGUMENTS
% data      Dim by nsubj, data which has nan where there is missing data.
% mask      A mask of the data which is made of 1s and 0s. 1s for where
%           there is data and nans for where there is no data. Default is
%           taken to be the mask with 1s everywhere.
%--------------------------------------------------------------------------
% OUTPUT
% fwhm_est      An estimate of the fwhm in each of the directions.
% Lambda_est    An estimate of the covariance matrix of the partial
%               derivatives.
% sigma_est     An estimate of the smoothness in terms of sigma.
%--------------------------------------------------------------------------
% EXAMPLES
% Dim = 160;
% nsubj = 20;
% noise = noisegen(Dim, nsubj, 6);
% est_smooth(noise)
%
% Dim = [250,250];
% nsubj = 1;
% noise = noisegen(Dim, nsubj, 20);
% est_smooth(noise)
% 
% Dim = [250,250];
% nsubj = 20;
% noise = noisegen(Dim, nsubj, 4);
% est_smooth(noise)
% 
% Dim = [250,250];
% nsubj = 100;
% noise = noisegen(Dim, nsubj, 6);
% est_smooth(noise)
% 
% Dim = [91,109,91];
% nsubj = 100;
% noise = noisegen(Dim, nsubj, 6);
% est_smooth(noise) %Gets around 6.14
% 
% % 3D example with a mask
% Dim = [91,109,91];
% nsubj = 20;
% noise = noisegen(Dim, nsubj, 3);
% mask = imgload('MNImask');
% est_smooth(noise, mask) %Gets around 3.26
%--------------------------------------------------------------------------
% AUTHORS: Samuel Davenport and Fabian Telschow
%--------------------------------------------------------------------------

%% %% Setup
%--------------------------------------------------------------------------
size_f = field.fieldsize;
D      = field.D;

if D >= 4
    error('Est_smooth only works in 1,2 and 3 dimensions')
end

nsubj = field.fibersize;

Y    = field.field;
mask = field.mask;

voxdim = get_dx( field );

%Make the zero-entries of the mask nan.
mask = zero2nan(mask);

% Index to remove cases for different dimensions
index = repmat( {':'}, 1, D );
nVox  = sum( ~isnan( mask(:) ) );

%% Standardize and Mask
if ~exist('df', 'var')
    Y = ( Y - mean( Y, length( size_f ) ) ) .* mask;
    df = 1;
else
    % Mask the data
    Y = Y .* mask;
end

% Estimate the variance across the image
var_est = sum( Y( ~isnan(Y) ).^2 ) / ( nVox * (nsubj - df) );

% Scale inut data by standard deviation
Y = Y / sqrt( var_est );

%% Estimate Lambda Matrix and FWHMs
Lambda_est      = zeros( D );
fwhm_est_kiebel = zeros( 1, D );

Xderivmate = diff( Y, 1, 1 ) / voxdim(1);
tmp        = ~isnan( Xderivmate( index{:}, 1 ) );
denom      = sum( tmp(:) ) * ( nsubj - df ); %The first half of this
% is the number of voxels in a given subject where Xderivmate is not nan.
%Since we have nsubj subjects we have to multiply this by nsubj - 1 to get
%the total number of voxels. Note its nsubj - 1 rather than nsubj since by
%subtracting the mean we have induced dependence across subjects.
%This is based off of Worsley's 1992 brain paper. Need to check how the
%derivation works to make things unbiased.

Lambda_est(1,1)    = sum( Xderivmate( ~isnan( Xderivmate ) ).^2 ) / denom;
fwhm_est_kiebel(1) = sqrt( 4 * log(2) / Lambda_est(1,1) );

if D > 1
    Yderivmate = diff( Y, 1, 2 ) / voxdim(2);
    tmp        = ~isnan( Yderivmate( index{:}, 1 ) );
    denom      = sum( tmp(:) ) * ( nsubj - df ); %The number of non-nans
    Lambda_est(2,2) = sum( Yderivmate( ~isnan( Yderivmate ) ).^2 ) / denom;
    fwhm_est_kiebel(2) = sqrt( 4 * log(2) / Lambda_est(2,2) );
end
if D > 2
    Zderivmate = diff( Y, 1, 3 ) / voxdim(3);
    tmp        = ~isnan( Zderivmate( index{:}, 1 ) );
    denom      = sum( tmp(:) ) * ( nsubj - df );  
    Lambda_est(3,3) = sum( Zderivmate( ~isnan(Zderivmate) ).^2 ) / denom;
    fwhm_est_kiebel(3) = sqrt( 4 * log(2) / Lambda_est(3,3) );
end

sigma_est = fwhm_est_kiebel / ( sqrt( 8 * log(2) ) );

sigma_est_forman = sqrt( -1 ./ 4 ./ log( 1 - diag(Lambda_est) / 2 ) );
fwhm_est_forman  = sigma2FWHM( sigma_est_forman );

end