function out = Gkerderiv2(x, FHWM)
sigma2 = FWHM2sigma(FHWM)^2;
out = (-1/sigma2 + x.^2/sigma2^2).*exp(-x.^2/(2*sigma2))/sqrt(2*pi*sigma2);
end
