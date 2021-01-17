%% %% A file which calculates the marginal p-value distributions for smoothing
%% at different times.
%% Gaussianize validity demonstration
FWHM = 3; resadd = 0; params = ConvFieldParams(FWHM, resadd, 0);
niter_vec = [10000,1000,1000,100];
nvox_vec = [100,1000,10000,100000];

global ncfloc
for I = 1:length(nvox_vec)
    niters = niter_vec(I);
    nvox = nvox_vec(I)
    for nsubj = [20,50,100]
        nsubj
        rng(mod(FWHM,5) + nsubj)
        orig_pval_store = [];
        G_smoothbefore_pval_store = [];
        G_smoothafter_pval_store = [];
        for J = 1:niters
            Y = wfield(nvox, nsubj);
            Y_G = Gaussianize(Y);
            Y_Gsmooth = convfield(Y_G, params);
            smoothY = convfield(Y, params);
            orig_tfield = mvtstat(smoothY.field);
            G_smoothY = Gaussianize(smoothY);
            G_tfield_smoothbefore = mvtstat(G_smoothY.field);
            G_tfield_smoothafter = mvtstat(Y_Gsmooth.field);
            orig_pval_store = [orig_pval_store, 1 - tcdf(orig_tfield, nsubj-1)'];
            G_smoothbefore_pval_store = [G_smoothbefore_pval_store, 1 - tcdf(G_tfield_smoothbefore, nsubj-1)'];
            G_smoothafter_pval_store = [G_smoothafter_pval_store, 1 - tcdf(G_tfield_smoothafter, nsubj-1)'];
        end
        clf
        alpha = 0.05;
        plot([0,alpha], [0,alpha])
        hold on
        [F, X] = ecdf(orig_pval_store);
        plot(X, F)
        xlim([0,alpha])
        ylim([0,alpha])
        [G_F_smoothbefore, G_X_smoothbefore] = ecdf(G_smoothbefore_pval_store);
        plot(G_X_smoothbefore,G_F_smoothbefore)
        [G_F_smoothafter, G_X_smoothafter] = ecdf(G_smoothafter_pval_store);
        plot(G_X_smoothafter,G_F_smoothafter, '--')
        legend('y = x', 'Original', 'Gaussianized smooth before',...
                            'Gaussianized smooth after', 'Location', 'SE')
        xlabel('x'); ylabel('F(x)'); title('Comparing the p-value distributia')
        export_fig([ncfloc, 'Figures/Simulations/PvalPlots/pvalsplot_nsubj_', num2str(nsubj), '_nvox_', num2str(nvox), '_FWHM_', num2str(FWHM),'.pdf'])
        clf
        h = histogram(G_smoothbefore_pval_store);
        h.BinWidth = 0.01;
        export_fig([ncfloc, 'Figures/Simulations/PvalPlots/SB_histogram_nsubj_', num2str(nsubj), '_nvox_', num2str(nvox), '_FWHM_', num2str(FWHM),'.pdf'])
        h = histogram(G_smoothafter_pval_store);
        h.BinWidth = 0.01;
        export_fig([ncfloc, 'Figures/Simulations/PvalPlots/SA_histogram_nsubj_', num2str(nsubj), '_nvox_', num2str(nvox), '_FWHM_', num2str(FWHM),'.pdf'])
    end
end

% Gaussianize validity T
field_type = 'T'; field_params = 3;
FWHM = 3; resadd = 0; params = ConvFieldParams(FWHM, resadd, 0);
niter_vec = [10000,1000,1000,100];
nvox_vec = [100,1000,10000,100000];

global ncfloc
for I = 1:length(nvox_vec)
    niters = niter_vec(I);
    nvox = nvox_vec(I)
    for nsubj = [20,50,100]
        nsubj
        rng(mod(FWHM,5) + nsubj)
        orig_pval_store = [];
        G_smoothbefore_pval_store = [];
        G_smoothafter_pval_store = [];
        for J = 1:niters
            Y = wfield(nvox, nsubj, field_type, field_params);
            Y_G = Gaussianize(Y);
            Y_Gsmooth = convfield(Y_G, params);
            smoothY = convfield(Y, params);
            orig_tfield = mvtstat(smoothY.field);
            G_smoothY = Gaussianize(smoothY);
            G_tfield_smoothbefore = mvtstat(G_smoothY.field);
            G_tfield_smoothafter = mvtstat(Y_Gsmooth.field);
            orig_pval_store = [orig_pval_store, 1 - tcdf(orig_tfield, nsubj-1)'];
            G_smoothbefore_pval_store = [G_smoothbefore_pval_store, 1 - tcdf(G_tfield_smoothbefore, nsubj-1)'];
            G_smoothafter_pval_store = [G_smoothafter_pval_store, 1 - tcdf(G_tfield_smoothafter, nsubj-1)'];
        end
        clf
        alpha = 0.05;
        plot([0,alpha], [0,alpha])
        hold on
        [F, X] = ecdf(orig_pval_store);
        plot(X, F)
        xlim([0,alpha])
        ylim([0,alpha])
        [G_F_smoothbefore, G_X_smoothbefore] = ecdf(G_smoothbefore_pval_store);
        plot(G_X_smoothbefore,G_F_smoothbefore)
        [G_F_smoothafter, G_X_smoothafter] = ecdf(G_smoothafter_pval_store);
        plot(G_X_smoothafter,G_F_smoothafter, '--')
        legend('y = x', 'Original', 'Gaussianized smooth before',...
                            'Gaussianized smooth after', 'Location', 'SE')
        xlabel('x'); ylabel('F(x)'); title('Comparing the p-value distributia')
        export_fig([ncfloc, 'Figures/Simulations/PvalPlots_T/pvalsplot_nsubj_', num2str(nsubj), '_nvox_', num2str(nvox), '_FWHM_', num2str(FWHM),'.pdf'])
        clf
        h = histogram(G_smoothbefore_pval_store);
        h.BinWidth = 0.01;
        export_fig([ncfloc, 'Figures/Simulations/PvalPlots_T/SB_histogram_nsubj_', num2str(nsubj), '_nvox_', num2str(nvox), '_FWHM_', num2str(FWHM),'.pdf'])
        h = histogram(G_smoothafter_pval_store);
        h.BinWidth = 0.01;
        export_fig([ncfloc, 'Figures/Simulations/PvalPlots_T/SA_histogram_nsubj_', num2str(nsubj), '_nvox_', num2str(nvox), '_FWHM_', num2str(FWHM),'.pdf'])
    end
end

%%
field_type = 'T'; field_params = 3;
niter_vec = [10000,1000,1000,100];
nvox_vec = [100,1000,10000,100000];

global ncfloc
for I = 1:length(nvox_vec)
    niters = niter_vec(I);
    nvox = nvox_vec(I)
    for nsubj = [20,50,100]
        nsubj
        rng(mod(FWHM,5) + nsubj)
        orig_pval_store = [];
        G_smoothbefore_pval_store = [];
        for I = 1:niters
            Y = wfield(nvox, nsubj, field_type, field_params);
            smoothY = convfield(Y, params);
            orig_tfield = mvtstat(smoothY.field);
            G_smoothY = Gaussianize(smoothY);
            G_tfield = mvtstat(G_smoothY.field);
            orig_pval_store = [orig_pval_store, 1 - tcdf(orig_tfield, nsubj-1)'];
            G_smoothbefore_pval_store = [G_smoothbefore_pval_store, 1 - tcdf(G_tfield, nsubj-1)'];
        end
        clf
        alpha = 0.05;
        plot([0,alpha], [0,alpha])
        hold on
        [F, X] = ecdf(orig_pval_store);
        plot(X, F)
        xlim([0,alpha])
        ylim([0,alpha])
        [G_F, G_X] = ecdf(G_smoothbefore_pval_store);
        plot(G_X,G_F)
        legend('y = x', 'Original p-values', 'Gaussianized p-values', 'Location', 'SE')
        xlabel('x'); ylabel('F(x)'); title('Comparing the p-value distributia')
        export_fig([ncfloc, 'Figures/Simulations/PvalPlots/T_pvalsplot_nsubj_', num2str(nsubj), '_nvox_', num2str(nvox), '_FWHM_', num2str(FWHM),'.pdf'])
        clf
        h = histogram(G_smoothbefore_pval_store);
        h.BinWidth = 0.01;
        export_fig([ncfloc, 'Figures/Simulations/PvalPlots/T_histogram_nsubj_', num2str(nsubj), '_nvox_', num2str(nvox), '_FWHM_', num2str(FWHM),'.pdf'])
        %     export_fig([ncfloc, 'Figures/Simulations/PvalPlots/nsubj_', num2str(nsubj)])
    end
end


%% Gaussianize validity demonstration
FWHM = 3; resadd = 0; params = ConvFieldParams(FWHM, resadd, 0);
niter_vec = [10000,1000,1000,100];
nvox_vec = [100,1000,10000,100000];

global ncfloc
for I = 1:length(nvox_vec)
    niters = niter_vec(I);
    nvox = nvox_vec(I)
    for nsubj = [20,50,100]
        nsubj
        rng(mod(FWHM,5) + nsubj)
        orig_pval_store = [];
        G_smoothbefore_pval_store = [];
        for J = 1:niters
            %         Y = wfield(nvox, nsubj, field_type, field_params);
            Y = wfield(nvox, nsubj);
            smoothY = convfield(Y, params);
            orig_tfield = mvtstat(smoothY.field);
            G_smoothY = Gaussianize(smoothY);
            G_tfield = mvtstat(G_smoothY.field);
            orig_pval_store = [orig_pval_store, 1 - tcdf(orig_tfield, nsubj-1)'];
            G_smoothbefore_pval_store = [G_smoothbefore_pval_store, 1 - tcdf(G_tfield, nsubj-1)'];
        end
        clf
        alpha = 0.05;
        plot([0,alpha], [0,alpha])
        hold on
        [F, X] = ecdf(orig_pval_store);
        plot(X, F)
        xlim([0,alpha])
        ylim([0,alpha])
        [G_F, G_X] = ecdf(G_smoothbefore_pval_store);
        plot(G_X,G_F)
        legend('y = x', 'Original p-values', 'Gaussianized p-values', 'Location', 'SE')
        xlabel('x'); ylabel('F(x)'); title('Comparing the p-value distributia')
        export_fig([ncfloc, 'Figures/Simulations/PvalPlots_SA/pvalsplot_nsubj_', num2str(nsubj), '_nvox_', num2str(nvox), '.pdf'])
        clf
        h = histogram(G_smoothbefore_pval_store);
        if nvox <= 1000
            h.BinWidth = 0.01;
        else
            h.BinWidth = 0.01;
        end
        export_fig([ncfloc, 'Figures/Simulations/PvalPlots_SA/histogram_nsubj_', num2str(nsubj), '_nvox_', num2str(nvox), '.pdf'])
        %     export_fig([ncfloc, 'Figures/Simulations/PvalPlots/nsubj_', num2str(nsubj)])
    end
end


field_type = 'T'; field_params = 3;
niter_vec = [10000,1000,1000,100];
nvox_vec = [100,1000,10000,100000];

global ncfloc
for I = 1:length(nvox_vec)
    niters = niter_vec(I);
    nvox = nvox_vec(I)
    for nsubj = [20,50,100]
        nsubj
        rng(mod(FWHM,5) + nsubj)
        orig_pval_store = [];
        G_smoothbefore_pval_store = [];
        for I = 1:niters
            %         Y = wfield(nvox, nsubj, field_type, field_params);
            Y = wfield(nvox, nsubj);
            smoothY = convfield(Y, params);
            orig_tfield = mvtstat(smoothY.field);
            G_smoothY = Gaussianize(smoothY);
            G_tfield = mvtstat(G_smoothY.field);
            orig_pval_store = [orig_pval_store, 1 - tcdf(orig_tfield, nsubj-1)'];
            G_smoothbefore_pval_store = [G_smoothbefore_pval_store, 1 - tcdf(G_tfield, nsubj-1)'];
        end
        clf
        alpha = 0.05;
        plot([0,alpha], [0,alpha])
        hold on
        [F, X] = ecdf(orig_pval_store);
        plot(X, F)
        xlim([0,alpha])
        ylim([0,alpha])
        [G_F, G_X] = ecdf(G_smoothbefore_pval_store);
        plot(G_X,G_F)
        legend('y = x', 'Original p-values', 'Gaussianized p-values', 'Location', 'SE')
        xlabel('x'); ylabel('F(x)'); title('Comparing the p-value distributia')
        export_fig([ncfloc, 'Figures/Simulations/PvalPlots_SA/pvalsplot_nsubj_', num2str(nsubj), '_nvox_', num2str(nvox), '.pdf'])
        clf
        h = histogram(G_smoothbefore_pval_store);
        if nvox <= 1000
            h.BinWidth = 0.01;
        else
            h.BinWidth = 0.01;
        end
        export_fig([ncfloc, 'Figures/Simulations/PvalPlots_SA/histogram_nsubj_', num2str(nsubj), '_nvox_', num2str(nvox), '.pdf'])
        %     export_fig([ncfloc, 'Figures/Simulations/PvalPlots/nsubj_', num2str(nsubj)])
    end
end


% %%
% % FWHM = 3; nsubj = 100;
% % load([where_davenpor, 'Data/Sim_Data/FWHM', num2str(FWHM), '_nsubj', num2str(nsubj)],...
% %                     'orig_pval_store', 'G_pval_store')
% clf
% alpha = 0.05;
% [F, X] = ecdf(orig_pval_store);
% plot(X, F)
% xlim([0,alpha])
% ylim([0,alpha])
% hold on
% [G_F, G_X] = ecdf(G_pval_store);
% plot(G_X,G_F)
% plot([0,alpha], [0,alpha])
% legend('Original p-values', 'Gaussianized p-values', 'y = x', 'Location', 'SE')
% xlabel('x'); ylabel('F(x)'); title('Comparing the p-value distributia')
% % global ncfloc
% export_fig([ncfloc,'Figures/Simulations/pvalue-plots/nsubj_', num2str(nsubj),...
%                             'FWHM', num2str(FWHM), '.pdf'], '-transparent')