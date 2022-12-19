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
% depth_text                       = {'0.0','0.5','1.0','1.5','2.0','2.5'};
% depth_val                        = [0.0,0.5,1.0,1.5,2.0,2.5];
Bins                             = 0:6;
Threshold                        = 0.3;
hems                             = {'lh','rh'};
sub_id                           = {'sub-01','sub-02','sub-03','sub-04'};
%ROIs                             = {'V1','V2','V3','V4'};
ROIs                             = {'V1','V2d','V2v','V3d','V3v','V4'};

% try
%     load([pRF_path 'results' filesep 'pRFsizes.mat'])
% catch
    PrfBinned                        = nan(length(sub_id),length(depth_text),length(Bins)-1,length(hems),length(ROIs));
    PrfEcc2                          = nan(length(sub_id),length(depth_text),length(hems),length(ROIs));
% end


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
                
                for roi = 1:length(ROIs)
                    try
                        Roi_file = [pRF_path ...
                            fullfile('data_FS',sub_id{sub},['fMRI/bar/depth_0.0/ROIs_tproc_' hems{hem} '_pRF_Gaussian' ...
                            filesep hems{hem} '_' ROIs{roi}])];
                        if ~exist([Roi_file '.label'])
                            error('ROIs not yet delineated!')
                        end
                    catch
                        Roi_file = '';
%                         if roi > 1
                            continue
%                         end
                    end
                    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                    disp([sub_id{sub} ', depth: ' depth_text{depth_idx} ', hemi: ' hems{hem} ', roi: ' ROIs{roi}])
                    [Res,~] = samsrf_plot(Srf, 'Sigma', SrfIv, 'Eccentricity', Bins, Roi_file, Threshold);

                    % bin results by eccentricity
                    BinnedRes = zeros(1,length(Bins)-1);
                    for ecc = 1:length(Bins)-1
                        BinnedRes(ecc) = mean(Res(Res(:,1) >= Bins(ecc) & Res(:,1) < Bins(ecc+1),2));
                    end
                    PrfBinned(sub,depth_idx,:,hem,roi) = BinnedRes;

                    % fit line, get estimate at 2deg
                    c = polyfit(Res(:,1),Res(:,2),1);
                    % Display evaluated equation y = m*x + b
                    disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
                    PrfEcc2(sub,depth_idx,hem,roi) = c(1)*2 + c(2);
                end
                
                %% convert Samsrf Srf file to FS label
                x0_idx    = 2;
                y0_idx    = 3;
                Sigma_idx = 4;
                R2_idx    = 8;
                Cmf_idx   = 10;
                if ~exist('FS_labels')
                    mkdir('FS_labels')
                end
                % x0
                if ~exist(['FS_labels' filesep hems{hem} '_x0.label'])
                    samsrf_srf2label(Srf, ['FS_labels' filesep hems{hem} '_x0'], x0_idx)
                end
                % y0
                if ~exist(['FS_labels' filesep hems{hem} '_y0.label'])
                    samsrf_srf2label(Srf, ['FS_labels' filesep hems{hem} '_y0'], y0_idx)
                end
                % pRF size (sigma)
                if ~exist(['FS_labels' filesep hems{hem} '_Sigma.label'])
                    samsrf_srf2label(Srf, ['FS_labels' filesep hems{hem} '_Sigma'], Sigma_idx)
                end
                % goodness of fit (R2)
                if ~exist(['FS_labels' filesep hems{hem} '_R2.label'])
                    samsrf_srf2label(Srf, ['FS_labels' filesep hems{hem} '_R2'], R2_idx)
                end
                % cortical magnification factor (CMF)
                if ~exist(['FS_labels' filesep hems{hem} '_Cmf.label'])
                    samsrf_srf2label(Srf, ['FS_labels' filesep hems{hem} '_Cmf'], Cmf_idx)
                end
                % polar angle
                Thrsh = 0.01;
                if ~exist(['FS_labels' filesep hems{hem} '_Polar.label'])
                    [Srf,Idx] = getPolarEccMaps(Srf,Thrsh,'Polar');
                    samsrf_srf2label(Srf, ['FS_labels' filesep hems{hem} '_Polar'], Idx)
                end
                % eccentricity
                Thrsh = [0.01 0 6.5];
                if ~exist(['FS_labels' filesep hems{hem} '_Eccentricity.label'])
                    [Srf,Idx] = getPolarEccMaps(Srf,Thrsh,'Eccentricity');
                    samsrf_srf2label(Srf, ['FS_labels' filesep hems{hem} '_Eccentricity'], Idx)
                end
                
                % copy labels to results dir
                results_dir = [pRF_path fullfile('results',sub_id{sub},['depth_' depth_text{depth_idx}])];
                if ~exist(results_dir)
                    mkdir(results_dir)
                end
                copyfile('FS_labels/*',results_dir)
                
                if depth_val(depth_idx)==0
                    copyfile(['ROIs_tproc_' hems{hem} '_pRF_Gaussian/*'],results_dir)
                end
            catch
            end
            
        end
    end
end
save([pRF_path 'results' filesep 'pRFsizes.mat'],'PrfBinned','PrfEcc2','depth_val','Bins','hems')

%% Plots
close all
% pRF size as a function of ecc
figure;
subs=[1 length(sub_id)];
sub_idx=1;
for sub = subs
    for hem = 1%:length(hems)
        subplot(1,2,sub_idx)
        for roi = 1:length(ROIs)
            y=squeeze(squeeze(PrfBinned(sub,1,:,hem,roi)));
            plot(Bins(2:end),y)
            hold on;
        end
    end
    sub_idx = sub_idx+1;
end
legend(ROIs)


% plot left/right hem separately
sub=4;
for hem = 1%:length(hems)
    figure(sub);
    if hem == 1
        title('Left hemisphere')
    else 
        title('Right hemisphere')
    end
    for ecc = 1:length(Bins)-1
        subplot(3,2,ecc)
        plot(depth_val,nanmean(squeeze(PrfBinned(:,:,ecc,hem,1)),1))
%         y=squeeze(squeeze(PrfBinned(sub,:,ecc,hem,1)));
%         plot(depth_val(~isnan(y)),y(~isnan(y)))
        title(['Ecc=' num2str(Bins(ecc)) '-' num2str(Bins(ecc)+1) 'deg'])
        xticks(depth_val)
        xticklabels(depth_val)
        xlabel('distance from GM/WM border (mm)')
        ylabel('pRF size (deg)')
    end
end


% plot central eccentricities
figure
roi = 1;
subs = 1;
depths=1:4;
for hem = 1%:length(hems)
    if hem == 1
        title('Left hemisphere')
    else 
        title('Right hemisphere')
    end
    
    if length(subs)>1
        plot(depth_val,nanmean(squeeze(nanmean(squeeze(squeeze(PrfBinned(subs,:,depths,hem,roi))),1)),2))
    else
        plot(depth_val,nanmean(squeeze(squeeze(PrfBinned(subs,:,depths,hem,roi))),2))
    end
    %title(['Ecc=' num2str(Bins(ecc)) '-' num2str(Bins(ecc)+1) 'deg'])
    xticks(depth_val)
    xticklabels(depth_val)
    xlabel('distance from GM/WM border (mm)')
    ylabel('pRF size (deg)')
    title([sub_id{subs} ', ROI: ' ROIs{roi}])
end

% plot both hems together
figure;
for ecc = 1:length(Bins)-1
    subplot(3,2,ecc)
    plot(depth_val,nanmean(squeeze(nanmean(squeeze(PrfBinned(:,:,ecc,:,1)),1)),2))
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


