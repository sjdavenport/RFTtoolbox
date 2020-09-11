%% Testing under the null
N = 20;
nvox = 10000;
field_type = 'L'; field_params = 3;
Y = wfield( nvox, N, field_type, field_params );
G_Y = Gaussianize(Y);

subplot(2,1,1)
histogram(1 - tcdf(mvtstat(Y.field), N-1), 'BinWidth', 0.05)
title('Pvalues - No Gaussianization')
subplot(2,1,2)
histogram(1 - tcdf(mvtstat(G_Y.field), N-1), 'BinWidth', 0.05)
title('Pvalues - Gaussianized')

%% Testing under the null with smoothing
N = 100;
FWHM = 3;
field_type = 'L'; 
nvox = 100000;
field_params = 3; % Only relevant if field_type is 'T' or 'L'
Y = wfield( nvox, N, field_type, field_params );
[smoothtstat, smooth_Y] = convfield_t(Y, FWHM);
G_Y = Gaussianize(Y);
G_Y_smoothbeforeG = Gaussianize(smooth_Y);
[~,G_smoothafterG] = convfield_t(G_Y, FWHM);

subplot(3,1,1)
histogram(1 - tcdf(smoothtstat.field, N-1), 'BinWidth', 0.05)
title('Pvalues - No Gaussianization with smoothing')
subplot(3,1,2)
histogram(1 - tcdf(mvtstat(G_Y_smoothbeforeG.field), N-1), 'BinWidth', 0.05)
title('Pvalues - Gaussianized - Smoothing Before')
subplot(3,1,3)
histogram(1 - tcdf(mvtstat(G_smoothafterG.field), N-1), 'BinWidth', 0.05)
title('Pvalues - Gaussianized - Smoothing After')
