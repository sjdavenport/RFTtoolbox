% Need to set rftboxloc to the right directory.
% export_fig is required to save the figures
rftboxloc = '/data/fireback/davenpor/davenpor/Toolboxes/RFTtoolbox/';

clf
pos_vector = [0,550,1000,600];
set(gcf, 'position', pos_vector)
set(0,'defaultAxesFontSize', 15);

noise = noisegen(160,20,6);
plot(mean(noise,2), 'linewidth', 2)
xlabel('Voxels')
ylabel('Mean over 20 subjects')
title('1D noise generation')

export_fig([rftboxloc, 'Figures/readme_1Dreal.png'], '-transparent')


clf
pos_vector = [0,550,1600,800];
set(0,'defaultAxesFontSize', 20);
set(gcf, 'position', pos_vector)

Dim = [100,100];
noise = noisegen(Dim, 20, 6);
noise_mean = mean(noise,3);
surf(noise_mean)
title('2D noise generation')
xlabel('x')
ylabel('y')

export_fig([rftboxloc, 'Figures/readme_2Dreal.png'], '-transparent')

clf
pos_vector = [0,550,1600,800];
set(0,'defaultAxesFontSize', 20);
set(gcf, 'position', pos_vector)

Sig = gensig([1.3,2], 3, [10,20], [100,150], {[40,30], [70,120]});
surf(Sig)
title('2D signal generation')
xlabel('x')
ylabel('y')


export_fig([rftboxloc, 'Figures/readme_signal.png'], '-transparent')
