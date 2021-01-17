%% %% A script that compares the null distributions

%% Gaussian
N = 50; nvox = 10000;
Y = wfield(nvox,N);
[G_Y, G_Y_nulldist] = Gaussianize(Y);
 
subplot(3,1,1)
histogram(G_Y_nulldist(:))
title('Null (X-muhat/sigmahat) distribution')
subplot(3,1,2)
histogram(Y.field(:))
title('Original data')
subplot(3,1,3)
histogram(G_Y.field(:))
title('Gaussianized data')

%% T
field_type = 'T'; 
field_params = 3; 
Y = wfield(nvox,N,field_type, field_params);
[G_Y, G_Y_nulldist, G_Y_nulldistnomeansub ] = Gaussianize(Y);

subplot(4,1,1)
histogram(G_Y_nulldist(:))
title('Null (X-muhat/sigmahat) distribution')
subplot(4,1,2)
histogram(G_Y_nulldistnomeansub(:))
title('(X/sigmahat) distribution')
subplot(4,1,3)
histogram(Y.field(:))
title('Original data')
subplot(4,1,4)
histogram(G_Y.field(:))
title('Gaussianized data')

%% Symmetric Skew
N = 100; nvox = 10000;
field_type = 's2'; 
Y = wfield(nvox,N,field_type);
[G_Y, G_Y_nulldist, G_Y_nulldistnomeansub ] = Gaussianize(Y);

subplot(2,1,1)
histogram(G_Y_nulldist(:))
title('Null (X-muhat/sigmahat) distribution')
subplot(2,1,2)
histogram(G_Y_nulldistnomeansub(:))
title('(X/sigmahat) distribution')

figure
subplot(2,1,1)
histogram(Y.field(:))
title('Original data')
subplot(2,1,2)
histogram(G_Y.field(:))
title('Gaussianized data')