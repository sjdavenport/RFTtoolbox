FWHM = 3;
sample_size = 50;

Dim = [5,5,5];
spfn = @(nsubj) normrnd(0,1,[Dim, nsubj]);

record_coverage( spfn, sample_size, FWHM, ones(Dim), 10000)

% 10000 iterations gives a coverage of 0.0713 for Dim = [5,5,5]!! and nsubj = 20. 
% Lattice one is 0.0518 though!!!! all are too high, indicates that
% nonstationarities are not being taken into account correctly. Especially
% sinec we're using a slightly smaller mask than we should be using.
%10000 iterations gives a coverage of 0.073 for Dim = [30,30,30]!!! Too high.
% (with the extra mask you get more like 0.08

%%
FWHM = 6;
sample_size = 10;

Dim = [5,5,5];
spfn = @(nsubj) normrnd(0,1,[Dim, nsubj]);

record_coverage( spfn, sample_size, FWHM, ones(Dim), 1000, 2)

%%
FWHM = 3;
sample_size = 20;

Dim = [50,50,50];
spfn = @(nsubj) normrnd(0,1,[Dim, nsubj]);

record_coverage( spfn, sample_size, FWHM, ones(Dim), 10000)

%10000 iterations gives a coverage of 0.073 for Dim = [30,30,30]!!! Too high.
% (with the extra mask you get more like 0.08


%%
nvox = 100;
xvals = 1:nvox;
plot(boot_smoothtfield_lat)
max_tfield_lat
max(tcf(top_lmlocs))

%%
FWHM = 3;
sample_size = 20;

Dim = [5,5];
spfn = @(nsubj) normrnd(0,1,[Dim, nsubj]);

mask = ones(Dim);

record_coverage( spfn, sample_size, FWHM, mask, 10000) 

%%
FWHM = 3;
sample_size = 20;

Dim = [30,30];
spfn = @(nsubj) normrnd(0,1,[Dim, nsubj]);

mask = ones(Dim);

record_coverage( spfn, sample_size, FWHM, mask, 1000) 
 