% Illusrating the edge effect
nvox = 100;
nsim = 1000;
smooth_fields = zeros(nsim, nvox);
FWHM = 10;
for I = 1:nsim
    I
   f = wfield(nvox);
   smoothf = convfield(f, FWHM);
   smooth_fields(I, :) = smoothf.field;
end

%%
dim = [20,20];
nsim = 1000;
smooth_fields = zeros([nsim, dim]);
FWHM = 3;
for I = 1:nsim
   I
   f = wfield(dim);
   smoothf = convfield(f, FWHM);
   smooth_fields(I, :, :) = smoothf.field;
end

%%
surf(squeeze(var(smooth_fields, 0,1)))

%%
FWHM = 10;
f = wfield(100);
smoothf = convfield(f, FWHM);
plot(smoothf)