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
k = 2.3; % Set k = 1 for a quadratic kernel and k = 2 for a quartic kernel
% k = 2.3 is very good because the 
FWHM = 4;
nvox = 50;

f = @(x) (1 - x.^2).^k;
fprime = @(x) 1;
fprime2 = @(x) 1;
% fprime = @(x) -2*k*x.*(1-x.^2).^(k-1);
% fprime2 = @(x) -2*k*(1-x.^2).^(k-1) + 4*k*(k-1)*x.^2.*(1-x.^2)^(k-2);

kern = SepKernel(1, f, fprime, fprime2, 1, 0);

lat_data = wfield([nvox,1]);
params_1 = ConvFieldParams(FWHM, 0);
smoothf = convfield(lat_data, params_1);

params = ConvFieldParams(kern, 15);
f = convfield(smoothf, params);
plot(smoothf, 'o-')
hold on
plot(f)


%% 1D linear (or close!) interpolation
clf
k = 1.1; % Set k = 1 for a quadratic kernel and k = 2 for a quartic kernel
nvox = 10; resadd = 101;
f = @(x) 1 - abs(x).^k;
fprime = @(x) 1;
fprime2 = @(x) 1;
% fprime = @(x) -2*k*x.*(1-x.^2).^(k-1);
% fprime2 = @(x) -2*k*(1-x.^2).^(k-1) + 4*k*(k-1)*x.^2.*(1-x.^2)^(k-2);

kern = SepKernel(1, f, fprime, fprime2, 1, 0);

lat_data = wfield([nvox,1]);
params_1 = ConvFieldParams(1, 0);
smoothf = convfield(lat_data, params_1);

params = ConvFieldParams(kern, resadd);
f = convfield(smoothf, params);
plot(smoothf, 'o-')
hold on
plot(f)

%% 1D sinc interpolation
clf
k = 1.1; % Set k = 1 for a quadratic kernel and k = 2 for a quartic kernel
nvox = 10; resadd = 11;
f = @(x) sinc(x);
fprime = @(x) 1;
fprime2 = @(x) 1;
% fprime = @(x) -2*k*x.*(1-x.^2).^(k-1);
% fprime2 = @(x) -2*k*(1-x.^2).^(k-1) + 4*k*(k-1)*x.^2.*(1-x.^2)^(k-2);

kern = SepKernel(1, f, fprime, fprime2, 1, 0);

lat_data = wfield([nvox,1]);
params_1 = ConvFieldParams(1, 0);
smoothf = convfield(lat_data, params_1);

params = ConvFieldParams(kern, resadd);
f = convfield(smoothf, params);
plot(smoothf, 'o-')
hold on
plot(f)

%% interp field - not quite sure what this does
clf
resadd = 4; nvox = 10; FWHM = 4;
lat_data = wfield([nvox,1]);
params_1 = ConvFieldParams(FWHM, 0);
smoothf = convfield(lat_data, params_1);

dx = 1/(1+resadd);
xvalues = 1:dx:length(smoothf.xvals{1});
y = interp(smoothf.field,resadd);
plot(smoothf, 'o-')
hold on
plot(f)



%% Compare the kernels
x = -1:0.001:1;
for k = 2:0.1:2.5
    f = @(x) (1 - x.^2).^k;
    plot(x, f(x));
    hold on
end
hold on
f = @(x) (1-abs(x));
plot(x, f(x));

%% Compare the kernels
x = -1:0.001:1;
powers = [1, 1.1, 1.25, 1.5, 1.75];
for I = 1:length(powers)
    f = @(x) (1-abs(x).^powers(I));
    plot(x, f(x));
    hold on
end

%% contruct a smooth interpolating kernel
ninv = 1/6;
x = [-ninv, ninv, 0];
y = [1 - ninv, 1 - ninv, 1];
xprime = [-ninv, ninv];
yprime = [1, -1];
coeffs = polyfit( x, y, xprime, yprime );

% test!
pp = @(x) piecewisepoly(x, [-1,-ninv,ninv,1], {[1,1], coeffs, [1,-1]});
x = -1:0.01:1;
plot(x, pp(x))

%% 1D close interpolation
clf
nvox = 10; resadd = 101;
FWHM = 4;
fprime = @(x) 1;
fprime2 = @(x) 1;

kern = SepKernel(1, pp, fprime, fprime2, 1, 0);

lat_data = wfield([nvox,1]);
params_1 = ConvFieldParams(FWHM, 0);
smoothf = convfield(lat_data, params_1);

params = ConvFieldParams(kern, resadd);
f = convfield(smoothf, params);
plot(smoothf, 'o-')
hold on
plot(f)

%% 1D tstat interp - this illustrates that the interpolated fields don't have constant
%% variance and if you divide by the standard deviation then the location of the 
%% peaks changes!!
clf
nvox = 10; resadd = 101;
nsubj = 20;
storefs = zeros(nsubj, 1021);
storefs_lat = zeros(nsubj, nvox);
FWHM = 4;
fprime = @(x) 1;
fprime2 = @(x) 1;

kern = SepKernel(1, pp, fprime, fprime2, 1, 0);

for I = 1:nsubj
    lat_data = wfield([nvox,1]);
    params_1 = ConvFieldParams(FWHM, 0);
    smoothf = convfield(lat_data, params_1);
    
    params = ConvFieldParams(kern, resadd);
    f = convfield(smoothf, params);
    
    storefs(I,:) = f.field;
    storefs_lat(I,:) = smoothf.field;
end

tf = mvtstat(storefs');
plot(mvtstat(storefs_lat'), '-o');
hold on
plot(f.xvals{1}, tf)
