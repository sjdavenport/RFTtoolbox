nvox = 100;
xvals_vecs = 1:nvox;
nsubj = 100;
lat_field = normrnd(0,1,nsubj,nvox);
FWHM = 5;
Kprime = @(x) GkerMVderiv(x,FWHM);
Kprime2 = @(x) GkerMVderiv2(x,FWHM);
truncation = 0;

peakloc = find_peak_locs(mean(lat_field), FWHM, 1:nvox);
% cov_val_deriv2 = cov(mate)
cov_mate = pointcov( peakloc, lat_field, FWHM )/nsubj;

mvnpdf([0,1]', [0,0]', cov_mate)

%%
u = -0.005
fun = @(f, fderiv2) abs(fderiv2).*mvnpdf([f,fderiv2], mu1', cov_mate);

mu1 = [0,0]';
Sigma1 = cov_mate;
fun2 = @(f, fderiv2) fderiv2*exp(-(1/2)*([f,fderiv2]'-mu1)'*inv(Sigma1)*([f,fderiv2]'-mu1))/sqrt(2*pi)^2/det(Sigma1)

fun(0,0.01)
fun2(0,0.01)

fun(u,0.01)

%%
% fun = @(x,y) x+y
integral2(fun, u, 1, -1, 0)

%%
MVint(fun, {[u,1],[-0.5,0]}, 0.0005)

%%
xvals_vecs = {-0.17:0.001:0.17,-0.04:0.001:0.04};
xvals = xvals2voxels(xvals_vecs);
z = fun(xvals(:,1), xvals(:,2));
zre = reshape(z, length(xvals_vecs{1}), length(xvals_vecs{2}));
surf(zre)

MVint(fun, {[-0.17,0.17], [-0.04, 0.04]}, 0.001)

%%

MVN = @(f, fderiv2) mvnpdf([f,fderiv2]', mu1, cov_mate);
MVN2 = @(f, fderiv2) exp(-(1/2)*([f,fderiv2]'-mu1)'*inv(cov_mate)*([f,fderiv2]'-mu1))/sqrt(2*pi)^2/det(cov_mate);

MVN(0,0)
MVN2(0,0)
%%
C = 100000;
x = [0,1]';
exp(-(1/2)*(x-mu1)'*inv(Sigma1)*(x-mu1))/sqrt(2*pi)/det(Sigma1)
