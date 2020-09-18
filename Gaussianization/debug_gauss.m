N = 5000;
xdim = 1000;
Mn = [zeros(1,xdim/2), linspace(-2,2,xdim/2)]/5; % Mean 
Sd = repmat([fliplr(linspace(0.1,2,xdim/4)),linspace(0.1,2,xdim/4)],1,2); % Standard deviation

%%
field_type = 's'; 
field_params = 3; % Only relevant if field_type is 'T' or 'L'
Y = Sd'.*(wfield(xdim,N, field_type, field_params).field);
% Y = Sd'.*(wfield(xdim,N, field_type, field_params).field) + Mn';
% Y = (wfield(xdim,N, field_type, field_params).field) - (wfield(xdim,N, field_type, field_params).field);
G_Y = Gaussianize(Y);

%%
field_type = 's'; 
field_params = 3; % Only relevant if field_type is 'T' or 'L'
% Y = Sd'.*(wfield(xdim,N, field_type, field_params).field) + Mn';
Y = Sd'.*((wfield(xdim,N, field_type, field_params).field) - (wfield(xdim,N, field_type, field_params).field)) + Mn';
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
a = wfield(1,N, field_type).field;
b = wfield(1,N, field_type).field;

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