clear all; clc; close all;
programs_path   = '/home/mayajas/Documents/programs';
scripts_path    ='/home/mayajas/Documents/project-00-7t-pipeline-dev/code/analysis-scripts/matlab';
pRF_path        = '/home/mayajas/scratch/project-00-7t-pipeline-dev/pRF/';

% add scripts path
addpath(scripts_path)

% add paths to required software
addpath(genpath(fullfile(programs_path,'samsrf')))

% define depth, sub, hemi variables
depth_text                       = {'-1.5','-1.0','-0.5','0.0','0.5','1.0','1.5','2.0','2.5'};
depth_val                        = [-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5];
Bins                             = 0:6;
Roi                              = '';
Threshold                        = 0.1;
hems                             = {'lh','rh'};
sub_id                           = {'sub-01','sub-02','sub-03','sub-04'};

try
    load([pRF_path 'results' filesep 'pRFsizes.mat'])
catch
    PrfBinned                        = nan(length(sub_id),length(depth_text),length(Bins)-1,length(hems));
    PrfEcc2                          = nan(length(sub_id),length(depth_text),length(hems));
end


for sub = 1:length(sub_id)
    for hem = 1:length(hems)
        for depth_idx = 1:length(depth_text)
            if ~isnan(PrfBinned(sub,depth_idx,:,hem))
                continue
            end
            
            cd([pRF_path fullfile('data_FS',sub_id{sub},['fMRI/bar/depth_' depth_text{depth_idx}])])
            try
                load(['postproc_' hems{hem} '_pRF_Gaussian.mat'])
                SrfIv = Srf;

                [Res,~] = samsrf_plot(Srf, 'Sigma', SrfIv, 'Eccentricity', Bins, Roi, Threshold);

                % bin results by eccentricity
                BinnedRes = zeros(1,length(Bins)-1);
                for ecc = 1:length(Bins)-1
                    BinnedRes(ecc) = mean(Res(Res(:,1) >= Bins(ecc) & Res(:,1) < Bins(ecc+1),2));
                end
                PrfBinned(sub,depth_idx,:,hem) = BinnedRes;

                % fit line, get estimate at 2deg
                c = polyfit(Res(:,1),Res(:,2),1);
                % Display evaluated equation y = m*x + b
                disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
                PrfEcc2(sub,depth_idx,hem) = c(1)*2 + c(2);

                %% convert Samsrf Srf file to FS label
                Sigma_idx = 4;
                R2_idx    = 8;
                Cmf_idx   = 10;
                if ~exist('FS_labels')
                    mkdir('FS_labels')
                end
                if ~exist(['FS_labels' filesep hems{hem} '_Sigma.label'])
                    samsrf_srf2label(Srf, ['FS_labels' filesep hems{hem} '_Sigma'], Sigma_idx)
                end
                if ~exist(['FS_labels' filesep hems{hem} '_R2.label'])
                    samsrf_srf2label(Srf, ['FS_labels' filesep hems{hem} '_R2'], R2_idx)
                end
                if ~exist(['FS_labels' filesep hems{hem} '_Cmf.label'])
                    samsrf_srf2label(Srf, ['FS_labels' filesep hems{hem} '_Cmf'], Cmf_idx)
                end

                results_dir = [pRF_path fullfile('results',sub_id{sub},['depth_' depth_text{depth_idx}])];
                if ~exist(results_dir)
                    mkdir(results_dir)
                end
                copyfile('FS_labels/*',results_dir)
            catch
            end
        end
    end
end
save([pRF_path 'results' filesep 'pRFsizes.mat'],'PrfBinned','PrfEcc2','depth_val','Bins','hems')

%% Plots
close all
% plot left/right hem separately
for hem = 1:length(hems)
    figure(hem);
    if hem == 1
        title('Left hemisphere')
    else 
        title('Right hemisphere')
    end
    for ecc = 1:length(Bins)-1
        subplot(3,2,ecc)
        plot(depth_val,nanmean(squeeze(PrfBinned(:,:,ecc,hem)),1))
        title(['Ecc=' num2str(Bins(ecc+1)) 'deg'])
        xticks(depth_val)
        xticklabels(depth_val)
        xlabel('distance from GM/WM border (mm)')
        ylabel('pRF size (deg)')
    end
end

% plot both hems together
figure;
for ecc = 1:length(Bins)-1
    subplot(3,2,ecc)
    plot(depth_val,nanmean(squeeze(nanmean(squeeze(PrfBinned(:,:,ecc,:)),1)),2))
    title(['Ecc=' num2str(Bins(ecc+1)) 'deg'])
    xticks(depth_val)
    xticklabels(depth_val)
    xlabel('distance from GM/WM border (mm)')
    ylabel('pRF size (deg)')
end

% plot pRF size at 2deg eccentricity only
figure;
subplot(3,1,1)
plot(depth_val,nanmean(PrfEcc2(:,:,1),1))
xlabel('distance from GM/WM border (mm)')
ylabel('pRF size (deg)')
title('Left hemisphere')
subplot(3,1,2)
plot(depth_val,nanmean(PrfEcc2(:,:,2),1))
xlabel('distance from GM/WM border (mm)')
ylabel('pRF size (deg)')
title('Right hemisphere')
subplot(3,1,3)
plot(depth_val,nanmean(squeeze(nanmean(PrfEcc2,1)),2))
xlabel('distance from GM/WM border (mm)')
ylabel('pRF size (deg)')
title('Both hemispheres')


