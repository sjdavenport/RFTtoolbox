N = 500;
xdim = 1000;

%%
field_type = 'L'; 
field_params = 3; % Only relevant if field_type is 'T' or 'L'
Y = (wfield(xdim,N, field_type, field_params).field);

% Y = Sd'.*(wfield(xdim,N, field_type, field_params).field-wfield(xdim,N, field_type, field_params).field) + Mn';

G_Y = Gaussianize(Y);
%%
G_Y2 = Gaussianize(G_Y);

plot(mvtstat(Y));hold on;plot(mvtstat(G_Y2.field));
legend('Original T-stat', 'Gaussianized T-stat', 'Location', 'Best')

%%
G_Y3 = Gaussianize(G_Y2);

plot(mvtstat(Y));hold on;plot(mvtstat(G_Y2.field));
legend('Original T-stat', 'Gaussianized T-stat', 'Location', 'Best')

%%
Mn = [zeros(1,xdim/2), linspace(-2,2,xdim/2)]/5; % Mean 
Sd = repmat([fliplr(linspace(0.1,2,xdim/4)),linspace(0.1,2,xdim/4)],1,2); % Standard deviation

clf
N = 100;
FWHM = 6;
field_type = 's2'; 
field_params = 3; % Only relevant if field_type is 'T' or 'L'
Y = Sd'.*(wfield(xdim,N, field_type, field_params).field) + Mn';
G_Y = Gaussianize(Y);
G_Y2 = Gaussianize(Gaussianize(G_Y));

subplot(3,1,1)
plot(mvtstat(Y));hold on;plot(mvtstat(G_Y.field));
legend('Original T-stat', 'Single Gauss', 'Location', 'NW')
subplot(3,1,2)
plot(mvtstat(Y));hold on;plot(mvtstat(G_Y.field), 'color', def_col('green'));
legend('Original T-stat', 'Double Gauss', 'Location', 'NW')
subplot(3,1,3)
plot(mvtstat(G_Y.field), 'color', def_col('red'));hold on;plot(mvtstat(G_Y2.field), 'color', def_col('green'));
legend('Single Gauss', 'Double Gauss', 'Location', 'NW')

