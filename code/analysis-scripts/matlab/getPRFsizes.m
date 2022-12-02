clear all; clc; close all;
programs_path   = '/home/mayajas/Documents/programs';
scripts_path    ='/home/mayajas/Documents/project-00-7t-pipeline-dev/scripts/';

% add scripts path
addpath(scripts_path)

% add paths to required software
addpath(genpath(fullfile(programs_path,'samsrf')))

try
    load('/home/mayajas/Documents/project-00-7t-pipeline-dev/results/pRFsizes.mat')
catch
    depth_text                       = {'m1.5','m1.0','m0.5','0.0','0.5','1.0','1.5','2.0','2.5'};
    depth_val                        = [-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5];
    Bins                             = 0:6;
    Roi                              = '';
    Threshold                        = 0.01;
    hems                             = {'lh','rh'};
    PrfBinned                        = zeros(length(depth_text),length(Bins)-1,length(hems));
    PrfEcc2                          = zeros(length(depth_text),length(hems));
    for hem = 1:length(hems)
        for depth_idx = 1:length(depth_text)
            cd(['/home/mayajas/scratch/project-00-7t-pipeline-dev/pRF/data_FS/sub-01/fMRI/bar/depth_' depth_text{depth_idx}])
            load(['postproc_' hems{hem} '_pRF_Gaussian.mat'])
            if depth_idx == 1
                SrfIv = Srf;
            end

            [Res,~] = samsrf_plot(Srf, 'Sigma', SrfIv, 'Eccentricity', Bins, Roi, Threshold);

            % bin results by eccentricity
            BinnedRes = zeros(1,length(Bins)-1);
            for ecc = 1:length(Bins)-1
                BinnedRes(ecc) = mean(Res(Res(:,1) >= Bins(ecc) & Res(:,1) < Bins(ecc+1),2));
            end
            PrfBinned(depth_idx,:,hem) = BinnedRes;

            % fit line, get estimate at 2deg
            c = polyfit(Res(:,1),Res(:,2),1);
            % Display evaluated equation y = m*x + b
            disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
            PrfEcc2(depth_idx,hem) = c(1)*2 + c(2);
        end
    end
    save('/home/mayajas/Documents/project-00-7t-pipeline-dev/results/pRFsizes.mat','PrfBinned','PrfEcc2','depth_val','Bins','hems')
end

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
        plot(depth_val,squeeze(PrfBinned(:,ecc,hem)))
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
    plot(depth_val,mean(squeeze(PrfBinned(:,ecc,:)),2))
    title(['Ecc=' num2str(Bins(ecc+1)) 'deg'])
    xticks(depth_val)
    xticklabels(depth_val)
    xlabel('distance from GM/WM border (mm)')
    ylabel('pRF size (deg)')
end

% plot pRF size at 2deg eccentricity only
figure;
subplot(3,1,1)
plot(depth_val,PrfEcc2(:,1))
xlabel('distance from GM/WM border (mm)')
ylabel('pRF size (deg)')
title('Left hemisphere')
subplot(3,1,2)
plot(depth_val,PrfEcc2(:,2))
xlabel('distance from GM/WM border (mm)')
ylabel('pRF size (deg)')
title('Right hemisphere')
subplot(3,1,3)
plot(depth_val,mean(PrfEcc2,2))
xlabel('distance from GM/WM border (mm)')
ylabel('pRF size (deg)')
title('Both hemispheres')


