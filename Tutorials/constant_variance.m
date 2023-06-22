FWHM = 5;
lat_data =  wfield(10000,10000);
smooth_data = convfield(lat_data, FWHM);

global PIloc
load([PIloc,'Variance/storevars'], 'allvars')

X = smooth_data.field;
% X = X./std(X,0,1);

%%
mean(std(X,0,1))

%%
mean(std(smooth_data.field/sqrt(allvars(FWHM)), 0, 1))