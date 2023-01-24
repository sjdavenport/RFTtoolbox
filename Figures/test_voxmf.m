clear all
close all

% rng(42)
mask = [20 20 20];
% Get lattice data
wn = wfield( mask, 1 );
% Input data
% Kernel: SepKernel object. Numeric creates a seperable isotropic Gaussian
% kernel with FWHM kernel
kernel_FWHM = 5;
% Derivative type: 0/1/2 supported giving the actual field or its first or
% second derivative
derivtype = 0; 
% resolution increase of the field, i.e. number of voxels added inbetween
% each voxel
resadd = 0;
% lat_mask: mask the lat_data before smoothing or not
lat_masked = true;
% Get the params object for a convfield
params = ConvFieldParams(kernel_FWHM*ones([1 3]),...
                         resadd,...
                         ceil( resadd / 2 ),...
                         false );
params2 = ConvFieldParams(kernel_FWHM*ones([1 3]),...
                         resadd + 1,...
                         ceil( (resadd + 1) / 2 ),...
                         false );

% Construct ConvField object.
cfield = convfield( wn, params, derivtype );
mask  = cfield.field > 0.08;
cfield.mask = mask;
cfield = Mask(cfield);

cfield2 = convfield( wn, params2, derivtype );
cfield2.mask = mask_highres( mask, resadd+1, ceil( (resadd + 1) / 2 ));

figure(1)
plot_voxmf( cfield )

figure(2)
plot_voxmf( cfield2 )