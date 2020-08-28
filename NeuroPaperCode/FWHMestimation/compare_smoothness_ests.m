function [forman, kiebel, conv] = compare_smoothness_ests( mask, nsubj, FWHM )
% compare_fwhm_ests( Dim, nsubj, FWHM ) takes a field of lattice data with 
% given mask and with nsubj subjects, smoothes it and then estimates Lambda
% and the FWHM using the different methods in order to compare them all.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  mask
%  nsubj
%  FWHM 
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% [forman, kiebel, conv] = compare_smoothness_ests( true([10,1]), 10, 0.5 )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'opt1', 'var' )
   % default option of opt1
   opt1 = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
pad = ceil(4*FWHM2sigma(FWHM));
mask = logical( pad_vals( mask, pad) );

lat_data = wnfield(mask, nsubj);

resadd = 0;
params = ConvFieldParams(FWHM, resadd);
params.lat_masked = false;
cfield = convfield(lat_data, params);
dcfield = convfield(lat_data, params, 1);

subset_index = cell(1, cfield.D);
for d = 1:cfield.D
    subset_index{d} = pad+1:size(mask,d)-pad;
end

Lambda = Riemmetric_est( cfield, dcfield );
conv_Lambda_est = mean(Lambda.field(subset_index{:}));
conv.fwhm = mean(sqrt(4*log(2)./Lambda.field(subset_index{:})));
conv.Lambda = conv_Lambda_est;

[forman.fwhm, kiebel.fwhm] = est_smooth(cfield.field(subset_index{:}, :));

kiebel.Lambda = 4*log(2)./kiebel.fwhm.^2;
forman.Lambda = 4*log(2)./forman.fwhm.^2;

end

