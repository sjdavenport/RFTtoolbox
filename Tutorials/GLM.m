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
error = wfield( dim, N);

% Signal
beta = constfield( [ 2, 2 ]', true( dim ) );

% Design matrix
x1    = ones( [ 1, N ] );
x2    = zeros( [ 1, N ] );
x2(1:4) = 1;

X = [ x1; x2 ]';

% Data
Y = X * beta + error;

% "Estimator"
hatX = (X'*X) \ X';

% Estimate
hatbeta = hatX * Y;

% residuals
residuals = Y - X * hatbeta;
df        = N - rank(X);

sigma2hat = sum(residuals.^2, 1) ./ df;
se = sqrt( (c * inv(X'* X) * c') .* sigma2hat);


% Compute t-field and residuals
T = c * hatbeta ./ se;
residuals = residuals ./ sigma2hat;