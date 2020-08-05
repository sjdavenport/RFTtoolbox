tot_nsubj = 700;

mask_sets = zeros([91,109, nsubj]);
RS_2D_data = zeros([91,109, nsubj]);

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

%%
D = 2; niters = 1000;
spfn = get_sample_fields( RS_2D_data, mask_sets, 2 );
FWHM = 3; sample_size = 20; resadd = 1;
params = ConvFieldParams( FWHM, resadd );
record_coverage( spfn, sample_size, params, niters)