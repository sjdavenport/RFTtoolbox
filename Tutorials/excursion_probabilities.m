nvox = 1000;
fwhm = 4;
dim = [nvox,1];

nsim = 1000;
c_0 = 1.3;
c_1 = 1.5;
c_0_count = 0;
c_1_count = 0;
lm_0 = 0;
lm_1 = 0;
for I = 1:nsim
    modul(I,100)
    lat_data = wfield(nvox);
    cfield = convfield(lat_data, fwhm);
    if max(cfield.field) > c_0
        c_0_count = c_0_count + 1;
    end
    if max(cfield.field) > c_1
        c_1_count = c_1_count + 1;
    end
    local_maxima = lmindices(cfield.field, 'all');
    field_vals = cfield.field(local_maxima);
    lm_0 = lm_0 + sum(field_vals > c_0);
    lm_1 = lm_1 + sum(field_vals > c_1);
end
c_0_count
c_1_count
lm_0/nsim
lm_1/nsim

%%
Lambda = FWHM2Lambda( fwhm, 1 );
expected_number_of_maxima_1 = (nvox*Lambda^(1/2))/(2*pi)*(exp(-c_1^2/2))
exp(-(c_1^2-c_0^2)/2)*lm_0/nsim
