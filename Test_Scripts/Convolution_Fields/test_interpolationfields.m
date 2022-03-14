%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests interpolation fields
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% 1D Interpolation kernels based on (1-x^2) and its powers
clf
k = 1; % Set k = 1 for a quadratic kernel and k = 2 for a quartic kernel
f = @(x) (1 - x.^2).^k;
fprime = @(x) 1;
fprime2 = @(x) 1;
% fprime = @(x) -2*k*x.*(1-x.^2).^(k-1);
% fprime2 = @(x) -2*k*(1-x.^2).^(k-1) + 4*k*(k-1)*x.^2.*(1-x.^2)^(k-2);

kern = SepKernel(1, f, fprime, fprime2, 1, 0);

nvox = 50;
lat_data = wfield([nvox,1]);
params_1 = ConvFieldParams(4, 0);
smoothf = convfield(lat_data, params_1);

params = ConvFieldParams(kern, 15);
f = convfield(smoothf, params);
plot(smoothf, 'o-')
hold on
plot(f)


%% 1D linear interpolation
clf
k = 1; % Set k = 1 for a quadratic kernel and k = 2 for a quartic kernel
f = @(x) 1 - abs(x);
fprime = @(x) 1;
fprime2 = @(x) 1;
% fprime = @(x) -2*k*x.*(1-x.^2).^(k-1);
% fprime2 = @(x) -2*k*(1-x.^2).^(k-1) + 4*k*(k-1)*x.^2.*(1-x.^2)^(k-2);

kern = SepKernel(1, f, fprime, fprime2, 1, 0);

nvox = 50;
lat_data = wfield([nvox,1]);
params_1 = ConvFieldParams(1, 0);
smoothf = convfield(lat_data, params_1);

params = ConvFieldParams(kern, 15);
f = convfield(smoothf, params);
plot(smoothf, 'o-')
hold on
plot(f)

%% Compare the kernels
x = -1:0.01:1;
for k = 1:5
    f = @(x) (1 - x.^2).^k;
    plot(x, f(x));
    hold on
end
hold on
plot(x, 1 - abs(x))

