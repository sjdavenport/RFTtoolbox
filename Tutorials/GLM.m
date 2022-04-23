% voxels to test
i = 1;
j = 1;

% Number of samples
N = 100;

% Contrast
c = [ 1, 0.5 ];

% Domain dimension
dim = [ 10 10 ];

% Grid values
x = ( 1:10 ) / 10;
[ xx, yy ] = meshgrid( x, x );

% Noise field
error_field   = wfield( dim, N);
errorpt = squeeze( error_field.field( i, j, : ) );

% Signal
beta = constfield( [ 2, 2 ]', true( dim ) );
betapt = squeeze( beta.field( i, j, : ) );

% Design matrix
x1    = ones( [ 1, N ] );
x2    = zeros( [ 1, N ] );
x2(1:4) = 1;

X = [ x1; x2 ]';

% Data
Y   = X * beta + error_field;
Ypt = X * betapt + errorpt;

% "Estimator"
hatX = (X'*X) \ X';

% Estimate
hatbeta = hatX * Y;
hatbetapt = hatX * Ypt;

sum( abs( squeeze( hatbeta.field(i,j,:)) - hatbetapt ) )

% residuals
residuals   = Y - X * hatbeta;
residualspt = Ypt - X * hatbetapt;

sum( abs( squeeze( residuals.field(i,j,:)) - residualspt ) )
df        = N - rank(X);

sigma2hat   = sum(residuals.^2, 1) ./ df;
sigma2hatpt = sum(residualspt.^2, 1) ./ df;
sum( abs( squeeze( residuals.field(i,j,:)) - residualspt ) )

% standard error of contrast estimate
se   = sqrt( (c * inv(X'* X) * c') .* sigma2hat );
sept = sqrt( (c * inv(X'* X) * c') .* sigma2hatpt );

% Compute t-field and standardized residuals
T = c * hatbeta ./ se;
residuals = residuals ./ sigma2hat;

Tpt = c * hatbetapt ./ sept;
residualspt = residualspt ./ sigma2hatpt;

sum( abs( squeeze( T.field(i,j,:)) - Tpt ) )
sum( abs( squeeze( residuals.field(i,j,:)) - residualspt ) )
