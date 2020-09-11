N = 500;
xdim = 1000;
Mn = [zeros(1,xdim/2), linspace(-2,2,xdim/2)]; % Mean 
Sd = repmat([fliplr(linspace(0.1,2,xdim/4)),linspace(0.1,2,xdim/4)],1,2); % Standard deviation

%%
field_type = 'T'; 
field_params = 3; % Only relevant if field_type is 'T' or 'L'
Y = Sd'.*(wfield(xdim,N, field_type, field_params).field) + Mn';
G_Y = Gaussianize(Y);

plot(mvtstat(Y));hold on;plot(mvtstat(G_Y.field));
legend('Original T-stat', 'Gaussianized T-stat', 'Location', 'Best')

%%
clf
N = 100;
FWHM = 10;
field_type = 'T'; 
field_params = 3; % Only relevant if field_type is 'T' or 'L'
Y = Sd'.*(wfield(xdim,N, field_type, field_params).field) + Mn';
[smoothtstat, smooth_Y] = convfield_t(Y, FWHM);
G_Y = Gaussianize(Y);
G_Y_smoothbeforeG = Gaussianize(smooth_Y);
[~,G_smoothafterG] = convfield_t(G_Y, FWHM);

subplot(3,1,1)
plot(smoothtstat);hold on;plot(mvtstat(G_Y_smoothbeforeG.field));
legend('Original T-stat', 'Gaussianized T-stat - Smoothed Before Gaussianization', 'Location', 'NW')
subplot(3,1,2)
plot(smoothtstat);hold on;plot(mvtstat(G_smoothafterG.field), 'color', def_col('green'));
legend('Original T-stat', 'Gaussianized T-stat - Smoothed After Gaussianization', 'Location', 'NW')
subplot(3,1,3)
plot(mvtstat(G_Y_smoothbeforeG.field), 'color', def_col('red'));hold on;plot(mvtstat(G_smoothafterG.field), 'color', def_col('green'));
legend('Gaussianized T-stat - Smoothed Before Gaussianization', 'Gaussianized T-stat - Smoothed After Gaussianization', 'Location', 'NW')