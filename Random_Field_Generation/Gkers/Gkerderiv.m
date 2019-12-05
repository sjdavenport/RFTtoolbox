function out = Gkerderiv(x, FHWM)
% val = GkerMV(1, 3)
% deriv = Gkerderiv(1, 3)
% h = 0.00001;
% valplushx = GkerMV(1+h, 3);
% (valplushx - val)/h

sigma2 = FWHM2sigma(FHWM)^2;
out = (-x/sigma2).*exp(-x.^2/(2*sigma2))/sqrt(2*pi*sigma2);
end
