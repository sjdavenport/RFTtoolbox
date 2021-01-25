%% Non-stationary comparisons of clustersize versus peak height
rng(101)
nvox = 10;
A = normrnd(0,1, [nvox,nvox]);
% A = A./diag(A)
Sigma = A*A';
Sigmavars = diag(Sigma);
global RFTboxloc
save([RFTboxloc, 'Random_Field_Generation/Nonstatnoisegeneration/nonstatvariables.mat'], 'Sigma', 'Sigmavars')

%%
% estimate variance and lambda (to make variance 1)
niters = 20000;
xvals = 1:0.001:10;
L = length(xvals);
field_store = zeros(niters,L);
resadd = 999;
enlarge = 0;
nvox = 10;

FWHM_vec = 3:7;
var_est = zeros(length(FWHM_vec), length(xvals));
for I = 1:length(FWHM_vec)
    FWHM = FWHM_vec(I);
    params = ConvFieldParams( FWHM, resadd, enlarge);
   
    for J = 1:niters
        modul(J,1000)
        lat_data = mvnrnd(zeros(1,nvox)', Sigma)./Sigmavars';
        field_store(J, :) = convfield(lat_data,params).field';
    end
    
    var_est(I,:) = var(field_store);
end

global RFTboxloc
save([RFTboxloc, 'Random_Field_Generation/Nonstatnoisegeneration/nonstatvar.mat'], 'var_est')

%%
global RFTboxloc
load([RFTboxloc, 'Random_Field_Generation/Nonstatnoisegeneration/nonstatvar.mat'], 'var_est')
load([RFTboxloc, 'Random_Field_Generation/Nonstatnoisegeneration/nonstatLambda.mat'], 'Lambda_est')

% Obtain Lambda
niters = 20000;
xvals = 1:0.001:10;
inter = 0.001;
L = length(xvals);
field_store = zeros(niters,L);
resadd = 999;
enlarge = 0;
nvox = 10;

FWHM_vec = 3:7;
Lambdascaled_est = zeros(length(FWHM_vec), length(xvals)-1);
for I = 1:length(FWHM_vec)
    FWHM = FWHM_vec(I);
    params = ConvFieldParams( FWHM, resadd, enlarge);
   
    for J = 1:niters
        modul(J,1000)
        lat_data = mvnrnd(zeros(1,nvox)', Sigma)./Sigmavars';
        field_store(J, :) = convfield(lat_data,params).field'./sqrt(var_est(I,:));
    end
    
    Lambdascaled_est(I,:) = var(diff(field_store,1,2)/inter);
end

global RFTboxloc
save([RFTboxloc, 'Random_Field_Generation/Nonstatnoisegeneration/nonstatLambda.mat'], 'Lambdascaled_est', 'Lambda_est')
