
load([data_path 'UKB_2D_randomized'])
if exist('im_store_2D', 'var')
    im_store = im_store_2D;
end

niters = 500;
sample_size = 100;
% nsubj=50, niters=5000
%    conv: 0.0428
%     lat: 0.0382
% finelat: 0.0418

% nsubj=200, niters=5000
%    conv: 0.0394
%     lat: 0.0348
% finelat: 0.0382

% randomized
% nsubj=100, niters=50
%    conv: 0.0360
%     lat: 0.0320
% finelat: 0.0360

figure
% nsubj=50, niters=5000
%    conv: 0.0508
%     lat: 0.0460
% finelat: 0.0496
% nsubj=200, niters=5000
%    conv: 0.0468
%     lat: 0.0422
% finelat: 0.0462

% randomized
% nsubj=100, niters=50
%    conv: 0.0480
%     lat: 0.0460
% finelat: 0.0480

figure
% nsubj=50, niters=5000
%    conv: 0.0528
%     lat: 0.0464
% finelat: 0.0514
% nsubj=200, niters=5000
%    conv: 0.0512
%     lat: 0.0452
% finelat: 0.0496
niters = 500;
sample_size = 100;
spfn_orig   = get_sample_fields( data2D_normalized, mask2D_eroded, D );
base_fields = spfn_orig(100);

data2df = base_fields.lat_data;
data2df.mask = mask2D_eroded;

[ L, L0, ~, cfields ] = LKC_latconv_est( data2df, params );

L_HPE = LKC_HP_est( cfields{1}, 1000, 1 );

% (This seems to work)
[ L; L_HPE.hatL' ]

% nsubj=50, niters=5000
%    conv: 0.0270
%     lat: 0.0192
% finelat: 0.0254
% nsubj=200, niters=5000
%    conv: 0.0282
%     lat: 0.0216
% finelat: 0.0264
% nsubj=50, niters=5000
%    conv: 0.0492
%     lat: 0.0378
% finelat: 0.0456
% nsubj=200, niters=5000
%    conv: 0.0446
%     lat: 0.0364
% finelat: 0.0426

% randomized
% nsubj=100, niters=500
%    conv: 0.0380
%     lat: 0.0280
% finelat: 0.0360

figure

% nsubj=50, niters=5000
%    conv: 0.0494
%     lat: 0.0282
% finelat: 0.0444
% nsubj=200, niters=5000
%    conv: 0.0520
%     lat: 0.0308
% finelat: 0.0454