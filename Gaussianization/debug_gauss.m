N = 5000;
xdim = 1000;
Mn = [zeros(1,xdim/2), linspace(-2,2,xdim/2)]/5; % Mean 
Sd = repmat([fliplr(linspace(0.1,2,xdim/4)),linspace(0.1,2,xdim/4)],1,2); % Standard deviation

%%
field_type = 's'; 
field_params = 3; % Only relevant if field_type is 'T' or 'L'
% @Sam: needed to separate, since my matlab version does not support
% object.field...., so your code depends on at least R2018b. We need to
% think about what we want to support!
f =  wfield( xdim, N, field_type, field_params );
Y = Sd' .* ( f.field );
% @Sam: note that this does the same on the fiber of the field class as
% what you did manually. Good practice would be if we keep variables in
% classes and only 
Y2 = Sd' .* f;
max(abs(Y(:)-Y2.field(:)))

% Y = Sd'.*(wfield(xdim,N, field_type, field_params).field) + Mn';
% Y = (wfield(xdim,N, field_type, field_params).field) - (wfield(xdim,N, field_type, field_params).field);
G_Y = Gaussianize(Y);

%%
field_type = 's'; 
field_params = 3; % Only relevant if field_type is 'T' or 'L'
% Y = Sd'.*(wfield(xdim,N, field_type, field_params).field) + Mn';
f1 = wfield( xdim, N, field_type, field_params );
f2 = wfield(xdim,N, field_type, field_params );
Y = Sd' .* ( f1.field - f2.field ) + Mn';
G_Y = Gaussianize(Y);

%%
a = (randn(500, 1).^2 - 1)/sqrt(2);
mean(a)
var(a)

%%
a = randn(500,1);
mean(a)
var(a)

%%
NGtstat = mvtstat(Y);
Gtstat = mvtstat(G_Y.field);

var(Gtstat(1:500))
var(NGtstat(1:500))

%%
N = 10000;
a = (randn(N, 1).^2 - 1)/sqrt(2);
Nb = 10000;
b = (randn(Nb, 1).^2 - 1)/sqrt(2);
mean(b)
c = b + 0.001;
vc = vec_ecdf( c, a );

subplot(2,1,1)
histogram(vc)

real_transform = zeros(1,length(c));
for I = 1:length(c)
    real_transform(I) = sum(a <= c(I));
end
real_transform = real_transform/length(real_transform);
subplot(2,1,2)
histogram(real_transform)

%%
N = 10000;
a = pearsrnd(0,1,1,4,N,1);
Nb = 10000;
b = pearsrnd(0,1,1,4,N,1);
mean(b)
c = b - 0.025;
vc = vec_ecdf( c, a );

subplot(2,1,1)
histogram(vc)

real_transform = zeros(1,length(c));
for I = 1:length(c)
    real_transform(I) = sum(a <= c(I));
end
real_transform = real_transform/length(real_transform);
subplot(2,1,2)
histogram(real_transform)

%%
N = 5000;
xdim = 1000;
field_type = 's2';
% Y = -(wfield(xdim,N, field_type).field);
f1 = wfield(1,N, field_type);
f2 = wfield(1,N, field_type);
a = f1.field;
b = f2.field;

c = b - 0.01;

vc = vec_ecdf( c, a );

subplot(2,1,1)
histogram(vc)

real_transform = zeros(1,length(c));
for I = 1:length(c)
    real_transform(I) = sum(a <= c(I));
end
real_transform = real_transform/length(real_transform);
subplot(2,1,2)
histogram(real_transform)

%%
N = 10000;
% a = (.^2 - 1)/sqrt(2);
b = (randn(N, 1).^2 - 1)/sqrt(2) - mean(b);
vc = vec_ecdf( a, b );
histogram(vc)