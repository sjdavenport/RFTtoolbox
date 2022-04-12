%--------------------------------------------------------------------------
%
%    Comparison of the LKC estimators
%
%--------------------------------------------------------------------------
%% prepare workspace
clear all
close all
%% 1D Isotropic field
%--------------------------------------------------------------------------
% Parameters
D    = 1;
FWHM = 6;
nsubj = 50;
T     = 100;
Msim = 100;
uvals = -4:0.5:4;
L = [1, NaN];

% Generate mask
pad    = ceil( 4*FWHM2sigma( FWHM ) );
mask = pad_vals( true( [ T, 1 ] ), pad, false );
lat_masked = false;

% Get high resolution theory
resThy = 301;
paramsThy = ConvFieldParams( FWHM, resThy, ceil(resThy/2), lat_masked );
LKC_thy_hr = LKC_wncfield_theory( mask, paramsThy );


%------------------
resaddVec = [0, 1, 3, 5];

L_voxM = cell([1, length(resaddVec)]);
L_conv = cell([1, length(resaddVec)]);
L_stat = cell([1, length(resaddVec)]);

L_HPE  = cell([1, length(resaddVec)]);
L_bHPE = cell([1, length(resaddVec)]);

L_reg_OLS  = cell([1, length(resaddVec)]);
L_breg_OLS = cell([1, length(resaddVec)]);
L_reg_SD   = cell([1, length(resaddVec)]);
L_breg_SD  = cell([1, length(resaddVec)]);
L_reg_SC   = cell([1, length(resaddVec)]);
L_breg_SC  = cell([1, length(resaddVec)]);
L_reg_PI   = cell([1, length(resaddVec)]);
L_breg_PI  = cell([1, length(resaddVec)]);

L_kiebel = cell([1, length(resaddVec)]);
L_forman = cell([1, length(resaddVec)]);

theoryLKC  = cell([1, length(resaddVec)]);


for r = 1:length(resaddVec)
    params = ConvFieldParams( FWHM, resaddVec(r), ceil(resaddVec(r)/2) , lat_masked );
    theoryLKC{r} = LKC_wncfield_theory( mask, params );

    lat_data =  wfield( mask, nsubj );

    % Convolution fields
    cfield  = convfield( lat_data, params );
    dcfield = convfield( lat_data, params, 1 );

    % Convolution fields estimators
    L_voxM{r}   = LKC_voxmfd_est( cfield, dcfield );
    L_conv{r} = LKC_latconv_est( lat_data, params );
    L_stat{r} = LKC_stationary_est( cfield, dcfield );

    % Convolution fields estimators
    L_HPE{r}      = LKC_HP_est( cfield, 1, 1 ).hatL;
    L_bHPE{r}     = LKC_HP_est( cfield, 200 ).hatL;
    L_reg_OLS{r}  = LKC_regress_est( cfield, uvals, "OLS", L, 1, 1 );
    L_breg_OLS{r} = LKC_regress_est( cfield, uvals, "OLS", L, 1, 200 );
    L_reg_SD{r}   = LKC_regress_est( cfield, uvals, "SD", L, 1, 1 );
    L_breg_SD{r}  = LKC_regress_est( cfield, uvals, "SD", L, 1, 200 );
    L_reg_SC{r}   = LKC_regress_est( cfield, uvals, "SC", L, 1, 1 );
    L_breg_SC{r}  = LKC_regress_est( cfield, uvals, "SC", L, 1, 200 );
    L_reg_PI{r}   = LKC_regress_est( cfield, uvals, "PI", L, 1, 1 );
    L_breg_PI{r}  = LKC_regress_est( cfield, uvals, "PI", L, 1, 200 );
    
    struct(...
        'theoryLKC_hr', LKC_thy_hr,...
        'theoryLKC_r', theoryLKC{r},...
        'L_voxM', L_voxM{r},...
        'L_conv', L_conv{r},...
        'L_stat', L_stat{r},...
        'L_HPE', L_HPE{r},...
        'L_bHPE', L_bHPE{r},...
        'L_reg_OLS', L_reg_OLS{r},...
        'L_breg_OLS', L_breg_OLS{r},...
        'L_reg_SD', L_reg_SD{r},...
        'L_breg_SD', L_breg_SD{r},...
        'L_reg_SC', L_reg_SC{r},...
        'L_breg_SC', L_breg_SC{r},...
        'L_reg_PI', L_reg_PI{r},...
        'L_breg_PI', L_breg_PI{r}...
        )
end


%% 2D Isotropic field
%--------------------------------------------------------------------------
% Parameters
D    = 2;
FWHM = repmat(6, [1 D]);
nsubj = 50;
T     = 30;
uvals = -4:0.5:4;
L = [1, NaN, NaN];

% Generate mask
pad    = ceil( 4*FWHM2sigma( FWHM ) );
mask = pad_vals( true( repmat(T, [1 D]) ), pad, false );
lat_masked = false;

% Get high resolution theory
resThy = 7;
paramsThy = ConvFieldParams( FWHM, resThy, ceil(resThy/2), lat_masked );
LKC_thy_hr = LKC_wncfield_theory( mask, paramsThy );


%------------------
resaddVec = [0, 1, 3, 5];

L_voxM = cell([1, length(resaddVec)]);
L_conv = cell([1, length(resaddVec)]);
L_stat = cell([1, length(resaddVec)]);

L_HPE  = cell([1, length(resaddVec)]);
L_bHPE = cell([1, length(resaddVec)]);

L_reg_OLS  = cell([1, length(resaddVec)]);
L_breg_OLS = cell([1, length(resaddVec)]);
L_reg_SD   = cell([1, length(resaddVec)]);
L_breg_SD  = cell([1, length(resaddVec)]);
L_reg_SC   = cell([1, length(resaddVec)]);
L_breg_SC  = cell([1, length(resaddVec)]);
L_reg_PI   = cell([1, length(resaddVec)]);
L_breg_PI  = cell([1, length(resaddVec)]);

L_kiebel = cell([1, length(resaddVec)]);
L_forman = cell([1, length(resaddVec)]);

theoryLKC  = cell([1, length(resaddVec)]);


for r = 1:length(resaddVec)
    params = ConvFieldParams( FWHM, resaddVec(r), ceil(resaddVec(r)/2) , lat_masked );
    theoryLKC{r} = LKC_wncfield_theory( mask, params );

    lat_data =  wfield( mask, nsubj+150 );

    % Convolution fields
    cfield  = convfield( lat_data, params );
    dcfield = convfield( lat_data, params, 1 );

    % Convolution fields estimators
    L_voxM{r} = LKC_voxmfd_est( cfield, dcfield );
    L_conv{r} = LKC_latconv_est( lat_data, params );
    L_stat{r} = LKC_stationary_est( cfield, dcfield );

    % Convolution fields estimators
    L_HPE{r}      = LKC_HP_est( cfield, 1, 1 ).hatL;
    L_bHPE{r}     = LKC_HP_est( cfield, 200 ).hatL;
    L_reg_OLS{r}  = LKC_regress_est( cfield, uvals, "OLS", L, 1, 1 );
    L_breg_OLS{r} = LKC_regress_est( cfield, uvals, "OLS", L, 1, 400 );
    L_reg_SD{r}   = LKC_regress_est( cfield, uvals, "SD", L, 1, 1 );
    L_breg_SD{r}  = LKC_regress_est( cfield, uvals, "SD", L, 1, 400 );
    L_reg_SC{r}   = LKC_regress_est( cfield, uvals, "SC", L, 1, 1 );
    L_breg_SC{r}  = LKC_regress_est( cfield, uvals, "SC", L, 1, 400 );
    L_reg_PI{r}   = LKC_regress_est( cfield, uvals, "PI", L, 1, 1 );
    L_breg_PI{r}  = LKC_regress_est( cfield, uvals, "PI", L, 1, 400 );
    
    struct(...
        'theoryLKC_hr', LKC_thy_hr,...
        'theoryLKC_r', theoryLKC{r},...
        'L_voxM', L_voxM{r},...
        'L_conv', L_conv{r},...
        'L_stat', L_stat{r},...
        'L_HPE', L_HPE{r}',...
        'L_bHPE', L_bHPE{r}',...
        'L_reg_OLS', L_reg_OLS{r},...
        'L_breg_OLS', L_breg_OLS{r},...
        'L_reg_SD', L_reg_SD{r},...
        'L_breg_SD', L_breg_SD{r},...
        'L_reg_SC', L_reg_SC{r},...
        'L_breg_SC', L_breg_SC{r},...
        'L_reg_PI', L_reg_PI{r},...
        'L_breg_PI', L_breg_PI{r}...
        )
end
