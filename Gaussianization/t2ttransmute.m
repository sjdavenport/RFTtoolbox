niters = 100000;
nsubj = 20;

meansubdist = zeros(1,niters);
targetdist = zeros(1,niters);
for I = 1:niters
    modul(I,1000)
    data1 = normrnd(0,1,1,nsubj);
    meansubdist(I) = (data1(1) - mean(data1))/std(data1);
    data2 = normrnd(0,1,1,nsubj);
    targetdist(I) = data2(1)/std(data2);
end
%%
h1 = histogram(meansubdist)
hold on
h2 = histogram(targetdist)
h2.BinWidth = h1.BinWidth

%%
niters = 100000;

for nsubj = 10:10:100
    meansubdist = zeros(1,nsubj*niters);
    targetdist = zeros(1,nsubj*niters);
    for I = 1:niters
        modul(I,1000)
        data1 = normrnd(0,1,1,nsubj);
        meansubdist((I-1)*nsubj+1:I*nsubj) = (data1 - mean(data1))/std(data1);
        data2 = normrnd(0,1,1,nsubj);
        targetdist((I-1)*nsubj+1:I*nsubj) = data2/std(data2);
    end
    save(['C:\Users\12SDa\davenpor\Data\RestingStateData\EklundData\Gaussianize_transmutes\GT2_nsubj_', ...
                                           num2str(nsubj)], 'meansubdist', 'targetdist')
end
%%
h1 = histogram(meansubdist)
hold on
h2 = histogram(targetdist)
h2.BinWidth = h1.BinWidth

%%
for nsubj = 10:10:100
    nsubj
    load(['C:\Users\12SDa\davenpor\Data\RestingStateData\EklundData\Gaussianize_transmutes\GT2_nsubj_', ...
        num2str(nsubj)], 'meansubdist', 'targetdist')
    vc_dist = vec_ecdf( targetdist, meansubdist );
    save(['C:\Users\12SDa\davenpor\Data\RestingStateData\EklundData\Gaussianize_transmutes\vc2_nsubj_',...
        num2str(nsubj)], 'vc_dist')
end

%%
nsubj = 50;
load(['C:\Users\12SDa\davenpor\Data\RestingStateData\EklundData\Gaussianize_transmutes\GT_nsubj_', ...
        num2str(nsubj)], 'meansubdist', 'targetdist')
h1 = histogram(meansubdist)
hold on
h2 = histogram(targetdist)
h2.BinWidth = h1.BinWidth

%%
nsubj = 30
load(['C:\Users\12SDa\davenpor\Data\RestingStateData\EklundData\Gaussianize_transmutes\vc_nsubj_',...
        num2str(nsubj)], 'vc_dist')
h = histogram(vc_dist)
h.BinWidth = 0.005
    
    