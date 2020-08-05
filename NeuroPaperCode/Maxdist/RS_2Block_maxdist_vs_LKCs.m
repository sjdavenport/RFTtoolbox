%%% Compare the maxima againt the EC curves obtained from the LKCs
FWHM_vec = 2:6;
nsubj_vec = [10,20,25,50,75];

global RFTboxloc
maxdistloc = [RFTboxloc, 'NeuroPaperCode/Maxdist/'];
for I = 1:length(FWHM_vec)
    FWHM = FWHM_vec(I);
    for J = 1:length(nsubj_vec)
       nsubj = nsubj_vec(I);
       load([maxdistloc, 'md_FWHM_', num2str(FWHM), 'nsubj_',num2str(nsubj)]);
    end
end

%%
FWHM = 6;
nsubj = 75;
MNImask = imgload('MNImask');
L0 = EulerChar(MNImask, 0.5, 3);

load([maxdistloc, 'md_FWHM_', num2str(FWHM), 'nsubj_',num2str(nsubj)]);
load([maxdistloc, 'LKC_estimates_MNImask'])

maxima = max(conv_maxima,[],1);
[ curve, x] = maxECcurve(maxima);
plot(x, curve)
hold on

LKCs = L_store_50(:,FWHM)';
LKC_est = [0, LKCs(2), LKCs(3)];
EEC_conv = EEC( x, L_store_50(:,FWHM)', L0, 'T', nsubj-1 )
plot(x,EEC_conv)
legend('maxima', 'EEC')