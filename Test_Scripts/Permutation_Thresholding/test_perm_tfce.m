%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the perm_tfce function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% Simple 1D example


%% %% 2D Examples
%% Simple 2D example
dim = [50,50]; nsubj = 50; FWHM = 0;
Sig = 0.25*peakgen(1, 10, 8, dim);
data = wfield(dim, nsubj).field + Sig;
threshold_tfce = perm_tfce(data, ones(dim));
tstat_orig = mvtstat(data);
tfce_tstat = tfce(tstat_orig);
threshold_orig = perm_tfce(data, ones(dim));

%%
dim = [50,50];
Sig = 0.35*square_signal(dim, 10, {[25,14], [25,36]} );
imagesc(Sig)

%% TFCE vs clustersize
connectivity_criterion = 8; H = 2; E = 0.5;
dim = [50,50]; nsubj = 50; 
% Sig = 0.5*peakgen(1, 10, 8, dim);
% Sig = zeros(dim); Sig(25:26,25) = 3;
% Sig = 0.35*square_signal(dim, 4, {[25,20], [25,30]} );
Sig = 0.35*square_signal(dim, 10, {[25,14], [25,36]} );
noise = wfield(dim, nsubj).field;
% FWHM = 0; 
FWHM = 6;
smooth_noise = fconv(noise, FWHM, 2);
smooth_noise = smooth_noise/sqrt(var(smooth_noise(:)));
data = smooth_noise + Sig;
threshold_tfce = perm_tfce(data, ones(dim), H, E, connectivity_criterion);
tstat_orig = mvtstat(data);
tfce_tstat = tfce(tstat_orig, H, E, connectivity_criterion);
%%
% CDT = 2.3;
CDT = 3.1;
threshold_cluster = perm_cluster(data, ones(dim), CDT, connectivity_criterion);

[number_of_clusters, occurences, sizes, index_locations] = numOfConComps(tstat_orig, CDT, connectivity_criterion);
surviving_cluster_im = cluster_im( dim, index_locations, threshold_cluster );

subplot(2,2,1)
imagesc(tstat_orig)
title('Original t-stat')
% view([-3, 15.58])
subplot(2,2,2)
imagesc(tfce_tstat)
% view([-3, 15.58])
title('TFCE statistic')
subplot(2,2,3)
imagesc(surviving_cluster_im)
title('Cluster extext inference: CDT = 3.1')
subplot(2,2,4)
imagesc(tfce_tstat > threshold_tfce)
% imagesc(tfce_tstat > 0)
title('TFCE: H = 2, E = 0.5, h_0 = 0')
% axis square
fullscreen
saveim('TFCEspill')

%%
save('./extfcedat', 'data', 'Sig') %in Testing in mp

%%
subjectids = filesindir('C:\Users\12SDa\davenpor\Data\HCP\HCPContrasts\WM\');

WM_imgs = zeros([91,109,91,80]);
for I = 1:length(subjectids)
    modul(I,10)
    subject_file_loc = ['C:\Users\12SDa\davenpor\Data\HCP\HCPContrasts\WM\',...
                                    subjectids{I}, '\WM\Level2\cope11.nii.gz'];
    WM_imgs(:,:,:,I) = imgload(subject_file_loc);
end

%%
slice = 35;
data = squeeze(WM_imgs(:,:,slice,:));
MNImask = imgload('MNImask');
mask = MNImask(:,:,slice);
connectivity_criterion = 8; H = 2; E = 0.5;
threshold_tfce = perm_tfce(data, mask, H, E, connectivity_criterion);
tstat_orig = mvtstat(data);
tfce_tstat = tfce(tstat_orig, H, E, connectivity_criterion);
CDT = 2.3;
threshold_cluster = perm_cluster(data, mask, CDT, connectivity_criterion);

[number_of_clusters, occurences, sizes, index_locations] = numOfConComps(tstat_orig, CDT, connectivity_criterion);
surviving_cluster_im = cluster_im( size(mask), index_locations, threshold_cluster );

%%
% imagesc(surviving_cluster_im)
slice = 35;
subplot(1,2,1)
overlay_brain([0, 0, slice], 10, {surviving_cluster_im}, 'red', 0.75)
title('Cluster size inference')
subplot(1,2,2)
overlay_brain([0, 0, slice], 10, {tfce_tstat> threshold_tfce}, 'red', 0.75)
title('TFCE: H = 2, E = 0.5')
fullscreen
% saveim('HCPslice')