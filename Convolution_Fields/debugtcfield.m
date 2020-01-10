L = 100; nsubj = 40;
data = lat_mean + sqrt(node_var).*normrnd(0,1, nsubj, L);
tcf = @(tval) tcfield( tval, data', xvalues_at_voxels, FWHM );
tcf2 = @(tval) tcfield2( tval, data, xvalues_at_voxels, FWHM );

xvals = 1:60;
plot(xvals, tcf(xvals) )
hold on
pause
plot(xvals, tcf2(xvals) )