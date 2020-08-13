data_path = 'C:\Users\12SDa\davenpor\Data\RestingStateData\EklundData\';
load([data_path 'UKB_2D'])
load([data_path, 'UKB_masks_2D'])
n_masks = size(mask_store, 3);

MNImask = imgload('MNImask');
MNImask_2D = MNImask(:,:,45);
new_mask = (1 - dilate_mask(mvtstat(im_store) < -3.6, 1)).*dilate_mask(MNImask_2D, -2);
imagesc(new_mask)

%%
slice = 60;
new_mask_1D = logical(squeeze(new_mask(slice,:)));
data_1D_nr = squeeze(im_store(slice,:,:));

%% View subject data
% data_1D = fconv(data_1D, 3, 1);
for I = 15:96
    subplot(2,1,1)
    plot(data_1D_nr(I,:))
    subplot(2,1,2)
    vox_data = data_1D_nr(I,:);
    BW = (max(vox_data) - min(vox_data))/30;
    histogram(vox_data, 'BinWidth', BW)
    pause
end

%% View smoothed subject data
smoothed_data_1D = fconv(data_1D, 3, 1);
for I = 15:96
    subplot(2,1,1)
    plot(smoothed_data_1D(I,new_mask_1D))
    subplot(2,1,2)
    vox_data = smoothed_data_1D(I,new_mask_1D);
    BW = (max(vox_data) - min(vox_data))/30;
    histogram(vox_data, 'BinWidth', BW)
    pause
end

%% Eklund distributia
RSDmask = imgload('RSDmask_Beijing');
data = loadCF(1:198, 'B1', 10);

[~,bounded_mask] = mask_bounds(RSDmask);
data_1D = squeeze(data(:,45,45,:));
mask_1D = bounded_mask(:,45,45)

%% View subject data
% data_1D = fconv(data_1D, 3, 1);
for I = find(mask_1D,1, 'first'):find(mask_1D,1, 'last')
    subplot(2,1,1)
    vox_data = data_1D(I,:);
    plot(vox_data)
    subplot(2,1,2)
    BW = (max(vox_data) - min(vox_data))/30;
    histogram(vox_data, 'BinWidth', BW)
    pause
end

%% Generate smoothed and averaged subject data
D = 1; spfn = get_sample_fields( data_1D, mask_1D, D );
nvox = size(data_1D,1); niters = 1000;
mean_store = zeros(niters, sum(mask_1D));
nsubj = 100; 
resadd = 0; FWHM = 3;
params = ConvFieldParams(FWHM, resadd);
for I = 1:niters
    modul(I)
    lat_data = spfn(nsubj).lat_data;
    smooth_data = convfield(lat_data, params);
    mean_store(I,:) = mean(lat_data.field(mask_1D,:), 2);
end


%% View smoothed and average data
for J = find(mask_1D,1, 'first'):find(mask_1D,1, 'last')
    clf
    subplot(2,1,1)
    plot(mean_store(:,J))
    subplot(2,1,2)
    p2v = mean_store(:,J)/sqrt(var(mean_store(:,J)));
    BW = 0.25;
    h = histogram(p2v, 'BinWidth', BW);
    x = h.BinLimits(1) + BW/2:BW:h.BinLimits(2);
    plot(x, h.Values/sum(h.Values))
    hold on
%     x = min(p2v):0.1:max(p2v);
    plot(x, normpdf(x, 0, 1)*BW);
    pause
end

%% 1D slice
D = 1; niters = 1000;
spfn = get_sample_fields( data_1D, mask_1D, D );
FWHM = 1; sample_size = 20; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
coverage1D = record_coverage( spfn, sample_size, params, niters)
% ECcurveanal( coverage1D, mask_slice, sample_size, 0.1 )
