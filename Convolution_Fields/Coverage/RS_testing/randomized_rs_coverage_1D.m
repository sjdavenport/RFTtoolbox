data_path = 'C:\Users\12SDa\davenpor\Data\RestingStateData\EklundData\';
load([data_path 'UKB_2D_randomized'])

MNImask = imgload('MNImask');
MNImask_2D = MNImask(:,:,45);

%% Histogram of the data 
slice = 60;
MNImask_1D = logical(squeeze(MNImask_2D(slice,:)));
data_1D = squeeze(im_store_2D(slice,:,:));
[tstat_1D,mu,std_1D] = mvtstat(data_1D);
% data_1D = fconv(data_1D, 10, 1);

%%
histogram(data_1D(logical(MNImask_1D),:), 'BinWidth', 0.25)
plot(mean(data_1D,2))

%% View subject data
% data_1D = fconv(data_1D, 3, 1);
for I = 15:96
    subplot(2,1,1)
    vox_data = data_1D(I,:);
    plot(vox_data)
    subplot(2,1,2)
    BW = (max(vox_data) - min(vox_data))/30;
    histogram(vox_data, 'BinWidth', BW)
    pause
end

%% View standardized data
p2v = mean_store(:,J)/sqrt(var(mean_store(:,J)));
BW = 0.25;
h = histogram(p2v, 'BinWidth', BW);
x = h.BinLimits(1) + BW/2:BW:h.BinLimits(2);
plot(x, h.Values/sum(h.Values))
hold on
%     x = min(p2v):0.1:max(p2v);
plot(x, normpdf(x, 0, 1)*BW);
pause

%% Plot tails
for I = 15:96
    clf
    vox_data = data_1D(I,:);
    std_data = vox_data/sqrt(var(vox_data));
    [a, x] = ecdf(abs(std_data));
    plot(x, log(1-a))
    hold on
    plot(x, log(1-normcdf(x)) + log(2))
    [b, bint] = regress(log(1-a(1:end-1)), x(1:end-1))
    hold on
    plot(x, b*x)
    pause
%     hold on
%     plot(x, -5*log(x))
%     pause
end

%%
y = normrnd(0,1,1,10000);
[a, x] = ecdf(abs(y));
plot(x, log(1-a))
hold on
plot(x, log(1-normcdf(x)) + log(2))
% plot(x, -log(x) - x.^2/2 -log(2), '-')

%% t
t_data = trnd(3, 1, 10000);
[a, x] = ecdf(abs(t_data));
plot(x, log(1-a))
hold on
plot(x, log(1-normcdf(x)) + log(2))
b = regress(log(1-a(1:end-1)), x(1:end-1));
hold on 
plot(x, b*x)
hold on
plot(x,-3*log(x))

%% Generate smoothed and averaged subject data
spfn = get_sample_fields( data_1D, MNImask_1D', D );
nvox = size(data_1D,1);
mean_store = zeros(niters, sum(MNImask_1D));
nsubj = 20;
resadd = 0; FWHM = 3;
params = ConvFieldParams(FWHM, resadd);
for I = 1:niters
    modul(I)
    lat_data = spfn(nsubj).lat_data;
    smooth_data = convfield(lat_data, params);
    mean_store(I,:) = mean(lat_data.field(MNImask_1D,:), 2);
end

%% View smoothed and average data
for J = 15:96
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

%% View smoothed subject data
smoothed_data_1D = fconv(data_1D, 3, 1);
for I = 15:96
    subplot(2,1,1)
    plot(smoothed_data_1D(I,:))
    subplot(2,1,2)
    histogram(smoothed_data_1D(I,:))
    pause
end

%% Plot smooth tails
FWHM = 10;
smoothed_data_1D = fconv(data_1D, FWHM, 1);
for I = 15:96
    clf
    vox_data = smoothed_data_1D(I,:);
    std_data = vox_data/sqrt(var(vox_data));
    [a, x] = ecdf(abs(std_data));
    plot(x, log(1-a))
    hold on
    plot(x, log(1-normcdf(x)) + log(2))
    [b, bint] = regress(log(1-a(1:end-1)), x(1:end-1))
    hold on
    plot(x, b*x)
    pause
%     hold on
%     plot(x, -5*log(x))
%     pause
end

%% tstat histogram
subplot(1,3,1)
histogram(tstat_1D(logical(MNImask_1D)), 'BinWidth', 0.25)
subplot(1,3,2)
plot(tstat_1D)
subplot(1,3,3)
plot(MNImask_1D)

%% Set Sample field function and parameters
D = 1; niters = 1000;
spfn = get_sample_fields( data_1D, MNImask_1D', D );

%% Record the coverage
sample_size = 200; resadd = 1; FWHM = 3;
params = ConvFieldParams ( repmat(FWHM,1,D), resadd );
coverage = record_coverage( spfn, sample_size, params, niters, 3)
ECcurveanal(coverage, MNImask_1D', sample_size, 0.1)

%% View fields
lat_data = spfn(75).lat_data;
FWHM = 20;
params = ConvFieldParams ( repmat(FWHM,1,D), resadd );
tcfield = convfield_t(lat_data, params);
subplot(1,2,1);
plot(tcfield)
subplot(1,2,2);
plot(MNImask_1D);

%% Record the coverage
D = 2; niters = 1000;
spfn = get_sample_fields( im_store_2D, MNImask_2D, D );

nsubj_vec = 20:20:110;
store_coverage = struct();
FWHM = 3; resadd = 11;
params = ConvFieldParams([FWHM, FWHM], resadd);

for I = 1:length(nsubj_vec)
    nsubj = nsubj_vec(I);
    store_coverage.(['nsubj_',num2str(nsubj)]) = record_coverage(spfn, nsubj, params, niters);
    save('./store_coverage_even', 'store_coverage')
end

%% EC curve analysis
ECcurveanal(coverage, MNImask_2D, sample_size, 0.1)
