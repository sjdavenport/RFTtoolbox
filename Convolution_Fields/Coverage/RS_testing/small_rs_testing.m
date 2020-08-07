tot_nsubj = 700;

mask_sets = zeros([91,109, tot_nsubj]);
RS_2D_data = zeros([91,109, tot_nsubj]);

global featruns
mask_dir = [featruns, 'RS_2Block_masks_warped/'];
directory = [featruns, 'RS_2Block_warped/'];

MNImask = imgload('MNImask');

for I = 1:tot_nsubj
    I
    mask = loadsubs(I, mask_dir, 0, ones([91,109,91]), 1);
    mask = (mask>0.5).*MNImask;
    im = loadsubs(I, directory, 0, ones([91,109,91]), 1);
    mask_sets(:,:,I) = mask(:,:,45);
    RS_2D_data(:,:,I) = im(:,:,45);
end

mask_sets = mask_sets > 0.5;

masks_1D = squeeze(mask_sets(:, 50, :));
RS_1D_data = squeeze(RS_2D_data(:, 50, :));

%%
D = 1; niters = 1000;
spfn = get_sample_fields( -RS_1D_data, masks_1D, D );
FWHM = 3; sample_size = 20; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
coverage = record_coverage( spfn, sample_size, params, niters)

%% 1D LKC calculations
lat_data = spfn(20).lat_data;
cfield  = convfield_Field( lat_data, FWHM, 0, resadd );
LKC_HP_est( cfield, 1, 1 )

resadd = 11;
lat_data = spfn(20).lat_data;
cfield  = convfield_Field( lat_data, FWHM, 0, resadd );
LKC_HP_est( cfield, 1, 1 )

cfield  = convfield_Field( lat_data, FWHM, 0, resadd );
dcfield = convfield_Field( lat_data, FWHM, 1, resadd );
[L_conv,L0_conv] = LKC_voxmfd_est( cfield, dcfield );

%%
D = 2; niters = 1000;
spfn = get_sample_fields( RS_2D_data, mask_sets, D );
FWHM = 3; sample_size = 20; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
record_coverage( spfn, sample_size, params, niters)