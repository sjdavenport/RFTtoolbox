clear all
close all

% General parameters
sigma = 5;
T     = 49;
FWHM  = sigma2FWHM( 5 );
nsubj = 20;
pad   = ceil(4*sigma);
dim   = [T T];
dimp  = dim + 2 * pad;
resadd  = 3;
enlarge = ceil( resadd / 2 );
xvals   = { (1:dimp(1)), (1:dimp(2)) };
% xvals   = { (1:dimp(1))/2, (1:dimp(2))/4 }; % different size of the box in the two directions 
theoryL = LKC_isogauss_theory( FWHM, [ T T ] )

% Mask
mask = true( dim );
mask = logical( pad_vals( mask, pad) );
%mask(:,2:4) = 0;


%% % Generate lattice data in simulation
%--------------------------------------------------------------------------
rng(1)
lat_data = wnfield( mask, nsubj);
% % Above lat_data field can equivalently obtained by
% rng(1)
% lat_data = wnfield( dim, nsubj) % mask is assumed to be all true. 
% lat_data.mask = mask

% Plot a slice at a point of a Field object with default
% (:, ceil( smask(2)/2 ), ..., ceil( smask(end)/2 ) )
figure(1), clf,
plot(lat_data)
figure(2), clf,
plot( lat_data, [ 12, NaN ] )

% Plot an image slice at a point of a Field object with default
% (:, :, ceil( smask(2)/2 ), ..., ceil( smask(end)/2 ) )
figure(3), clf
imagesc( lat_data )

% Different slice and subj
figure(4), clf
imagesc( lat_data, [ NaN, NaN ], 3 ) % @Fabian, this line doesn't work

% Want to change voxelsize?
lat_data.xvals = xvals;
% %Shorter in that case is
% rng(1)
% lat_data = wnfield( mask, nsubj, xvals );

%% % Generate lattice data in practice given obseved_data and a mask
%--------------------------------------------------------------------------
observed_data = lat_data.field;
% Create Field object from the observed array data and the mask
lat_data_obs = Field( observed_data, mask, xvals );
% By construction lat_data_obs is the same as lat_data. xvals is optional,
% if not provided voxels are assumed to have dimension 1 in all directions.

%% % Getting convolution fields
%--------------------------------------------------------------------------
% logical indicating whether the lat_data gets masked before smoothing or
% not,
lat_masked = true;

% Generate convolution fields from lattice data
cfield  = convfield_Field( lat_data, FWHM, 0, resadd, lat_masked, enlarge );
dcfield = convfield_Field( lat_data, FWHM, 1, resadd, lat_masked, enlarge );

figure, clf
imagesc( cfield )

figure, clf  
imagesc( Mask( cfield) )

% Mask lat_data manually for the old convfield
if lat_masked
    lat_data = Mask( lat_data );
end

[ smooth_data ] = convfield( lat_data.field, FWHM, resadd, 2,...
                                                         0, enlarge );
[ dsmooth_data ] = convfield( lat_data.field, FWHM, resadd, 2,...
                                                         1, enlarge );
                                                                            
% new convfield_Field function gives the same output as convfield, if xvals
% assumes all voxels have size 1. The other case was not supported before.
[ maxdiff( smooth_data, cfield.field ),...                                                    
  maxdiff( dsmooth_data, dcfield.field ) ]

%% % Getting the LKCs
%--------------------------------------------------------------------------
%%% Getting the LKCs
LKC_voxmfd_est( cfield, dcfield )

%%% The above function basically calls the following. You see the results
%%% are identical
% Estimate Riemannian metric induced by a convolution field.
Lambda = Lambda_est( cfield, dcfield );

% Obtain Voxmanifold object from the convfields.
% The property voxmfd.g = Lambda_est( cfield, dcfield ).
voxmfd = ConvField2VoxManifold( cfield, dcfield );

LKC_est( voxmfd ) 