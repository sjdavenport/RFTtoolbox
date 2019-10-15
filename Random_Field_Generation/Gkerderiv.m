function out = Gkerderiv(x, FHWM)
sigma2 = FWHM2sigma(FHWM)^2;
out = (-x/sigma2).*exp(-x.^2/(2*sigma2))/sqrt(2*pi*sigma2);
end
