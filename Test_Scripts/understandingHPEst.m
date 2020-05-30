% (A file for me to help under the HPE estimation)
dim   = [25 25];
nsubj = 30;
D     = length( dim );
mask  = true(dim);

% generate random noise
Y = noisegen( dim, nsubj, 5*2*sqrt(2*log(2)), 0 );
% compute HPE and bHPE
HPE  = LKCestim_HPE( Y, D, mask, 1, "C" );
bHPE = LKCestim_HPE( Y, D, mask, 5e3, "C" );

%%
HPE.hat1

mean(HPE.hat1, 2)
HPE.hatn

% (A file for me to help under the HPE estimation)
dim   = [10 10 10];
nsubj = 10;
D     = length( dim );
mask  = true(dim);

% generate random noise
Y = noisegen( dim, nsubj, 5*2*sqrt(2*log(2)), 0 );
% compute HPE and bHPE
HPE  = LKCestim_HPE( Y, D, mask, 1, "C" );

% In 3D, HPE.hatn is [L1hat, L2hat, L3hat], and L0 is 0.

%%
HPE.hat1

mean(HPE.hat1, 2)
HPE.hatn

cov(HPE.hat1')
HPE.Sigma_hat

%% EC evaluation
sY = size(Y);
ECs = EulerCharCrit( Y, D, mask );
ECs{1}
ECs{2}
%The first column gives the levels at which the EC of each curve changes
%value. The 2nd gives the EC at and above those levels. First row is always
%-Inf and 1, since everything is above -Inf so the EC becomes the EC of the
%containing set which is 1 here (as on a rectanguloid). The last row is
%always Infinifty

%% 1D EC evaluation
% (A file for me to help under the HPE estimation)
dim   = [1, 100];
nsubj = 30;
D     = length( dim );
mask  = true(dim);

% generate random noise
Y = noisegen( dim, nsubj, 5*2*sqrt(2*log(2)), 0 );
ECs = EulerCharCrit( Y, D, mask );
%%
plot(Y(:,:,1))
EC_change = zeros(size(ECs{1}(2:end-1, :)));
EC_change(:,1) = ECs{1}(2:end-1, 1);
difference = ECs{1}(2:end,2) - ECs{1}(1:end-1,2);
EC_change(:,2) = difference(2:end)

%% 2D EC evaluation
% (A file for me to help under the HPE estimation)
dim   = [25, 25];
nsubj = 30;
D     = length( dim );
mask  = true(dim);

% generate random noise
Y = noisegen( dim, nsubj, 5*2*sqrt(2*log(2)), 0 );
ECs = EulerCharCrit( Y, D, mask );
EC_change = zeros(size(ECs{1}(2:end-1, :)));
EC_change(:,1) = ECs{1}(2:end-1, 1);
difference = ECs{1}(2:end,2) - ECs{1}(1:end-1,2);
EC_change(:,2) = difference(2:end)
%%
firstimage = Y(:,:,1);
min(firstimage(:))
surf(Y(:,:,1))
%% 3D EC evaluation
% (A file for me to help under the HPE estimation)
dim   = [25, 25, 25];
nsubj = 10;
D     = length( dim );
mask  = true(dim);

% generate random noise
Y = noisegen( dim, nsubj, 5*2*sqrt(2*log(2)), 0 );
ECs = EulerCharCrit( Y, D, mask );
EC_change = zeros(size(ECs{1}(2:end-1, :)));
EC_change(:,1) = ECs{1}(2:end-1, 1);
difference = ECs{1}(2:end,2) - ECs{1}(1:end-1,2);
EC_change(:,2) = difference(2:end)
%%
firstimage = Y(:,:,1);
min(firstimage(:))
surf(Y(:,:,1))

%% 1D EC evaluation that doesn't work
% (A file for me to help under the HPE estimation)
dim   = 100;
nsubj = 30;
D     = length( dim );
mask  = true(dim);

% generate random noise
Y = noisegen( dim, nsubj, 5*2*sqrt(2*log(2)), 0 );
ECs = EulerCharCrit( Y, D, mask );