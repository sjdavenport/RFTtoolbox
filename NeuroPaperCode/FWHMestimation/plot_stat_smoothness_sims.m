%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script plots the stationary smoothness simulations
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D 
%% Load the data
global RFTboxloc
load([RFTboxloc, 'NeuroPaperCode/FWHMestimation/oneDstatsims_nsubj100']);
sim_types = {'forman', 'kiebel', 'conv'};
sim_names = {'Forman', 'Kiebel', 'Convolution'};

%% FWHM bias plot
startat = 1.5; start_index = find(FWHM_vec == startat);
plot(FWHM_vec(start_index:end), mean(forman.fwhm_ests(start_index:end,:),2) - FWHM_vec(start_index:end)')
hold on
plot(FWHM_vec(start_index:end), mean(kiebel.fwhm_ests(start_index:end,:),2) - FWHM_vec(start_index:end)')
plot(FWHM_vec(start_index:end), mean(conv.fwhm_ests(start_index:end,:),2) - FWHM_vec(start_index:end)')
legend(sim_names{1}, sim_names{2}, sim_names{3}, 'Location', 'NW')

%% FWHM var plot
startat = 1.5; start_index = find(FWHM_vec == startat);
plot(FWHM_vec(start_index:end), std(forman.fwhm_ests(start_index:end,:),0,2))
hold on
plot(FWHM_vec(start_index:end), std(kiebel.fwhm_ests(start_index:end,:),0,2))
plot(FWHM_vec(start_index:end), std(conv.fwhm_ests(start_index:end,:),0,2))
legend(sim_names{1}, sim_names{2}, sim_names{3}, 'Location', 'NW')

%% Lambda plot
startat = 1.5; start_index = find(FWHM_vec == startat);
plot(FWHM_vec(start_index:end), mean(forman.Lambda_ests(start_index:end,:),2))
hold on
plot(FWHM_vec(start_index:end), mean(kiebel.Lambda_ests(start_index:end,:),2))
plot(FWHM_vec(start_index:end), mean(conv.Lambda_ests(start_index:end,:),2))
legend(sim_names{1}, sim_names{2}, sim_names{3}, 'Location', 'NW')

%% Lambda var plot
startat = 1.5; start_index = find(FWHM_vec == startat);
plot(FWHM_vec(start_index:end), std(forman.Lambda_ests(start_index:end,:),0,2))
hold on
plot(FWHM_vec(start_index:end), std(kiebel.Lambda_ests(start_index:end,:),0,2))
plot(FWHM_vec(start_index:end), std(conv.Lambda_ests(start_index:end,:),0,2))
legend(sim_names{1}, sim_names{2}, sim_names{3}, 'Location', 'NW')