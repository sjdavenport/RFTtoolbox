%% 1D note that using nsubj = 10000 and high spacing we see that the true LKC is around 10.24
clc
nsubj = 1000;
lat_data = normrnd(0,1,nvox,nsubj);
mask = ones(1,nvox);

FWHM = 3;
D = 1;
spacing = 1;
LDcalc( lat_data, FWHM, D, mask, spacing)

% Compare to HPE/LKCGaussconv
Gker_param = FWHM2sigma(FWHM);
resAdd = floor(1/spacing);
pad = 0;
mask = ones(1,nvox);

% Convolution estimate
LKC_GaussConv( lat_data, Gker_param, D, resAdd, pad )
% LKCs_conv = LKCestim_GaussConv( Y, FWHM, D, resAdd, pad )

% HPE (remember need to smooth the fields before input for HPE!)
[~,smoothed_fields] = smoothtstat(lat_data, FWHM);

HPE  = LKCestim_HPE( smoothed_fields, 1, mask, 1);
HPE.hatn

%% 1D comparison
nvox = 2;
niters = 1000;
Gker_param = FWHM2sigma(FWHM);
pad = 0;
FWHM = 3;
D = 1;
spacing = 1;
resAdd = floor(1/spacing) -1;

nsubj = 50;

LDcalc_LKCs = zeros(1, niters);
LKCconv_LKCs = zeros(1, niters);
for I = 1:niters
    I
    lat_data = normrnd(0,1,nvox,nsubj);
    mask = ones(1,nvox);
    
    LDcalc_LKCs(I) = LDcalc( lat_data, FWHM, D, mask, spacing );
   
    % Convolution estimate
    LKCconv_LKCs(I) = LKC_GaussConv( lat_data, Gker_param, D, resAdd, pad );
end

LD_mean_LKC = mean(LDcalc_LKCs)
LD_var = var(LDcalc_LKCs)

LKCconv_mean_LKC = mean(LKCconv_LKCs)
LKCconv_var = var(LKCconv_LKCs)
