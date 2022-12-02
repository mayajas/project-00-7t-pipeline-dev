function getPRFsizes_curta()
clear all; clc; %close all;

project_name    ='project-00-7t-pipeline-dev';
programs_path   = '/home/mayaaj90/programs';
scripts_path    = ['/home/mayaaj90/scripts/' project_name filesep];
data_path       = ['/scratch/mayaaj90/' project_name filesep 'pRF/data_FS/sub-01/fMRI/bar/'];
res_path        = ['/scratch/mayaaj90/' project_name 'pRF/results'];

% add scripts path
addpath(scripts_path)

% add paths to required software
addpath(genpath(fullfile(programs_path,'samsrf')))

try
    load([res_path 'pRFsizes.mat'])
catch
    depth_text                       = {'m0.5','0.0','0.5','1.0','1.5','2.0','2.5'};
    depth_val                        = [-0.5,0.0,0.5,1.0,1.5,2.0,2.5];
    Bins                             = 0:6;
    Roi                              = '';
    Threshold                        = 0.3;
    Mode                             = 'Mean';
    Color                            = 'k';
    BootParams                       = [1000, 2.5, 97.5];
    hems                             = {'lh','rh'};
    Prfs                             = cell(length(hems),length(depth_text));
%     PrfBinned                        = zeros(length(depth_text),length(Bins)-1,length(hems));
%     PrfEcc2                          = zeros(length(depth_text),length(hems));
    for hem = 1:length(hems)
        for depth_idx = 1:length(depth_text)
            cd([data_path 'depth_' depth_text{depth_idx}])
            load(['postproc_' hems{hem} '_pRF_Gaussian.mat'])
            if depth_idx == 1
                SrfIv = Srf;
            end

            [Res,~] = samsrf_plot(Srf, 'Sigma', SrfIv, 'Eccentricity', Bins, Roi, Threshold, Mode, Color, BootParams);

            Prfs{hem,depth_idx} = Res;
%             % bin results by eccentricity
%             BinnedRes = zeros(1,length(Bins)-1);
%             for ecc = 1:length(Bins)-1
%                 BinnedRes(ecc) = mean(Res(Res(:,1) >= Bins(ecc) & Res(:,1) < Bins(ecc+1),2));
%             end
%             PrfBinned(depth_idx,:,hem) = BinnedRes;
% 
%             % fit line, get estimate at 2deg
%             c = polyfit(Res(:,1),Res(:,2),1);
%             % Display evaluated equation y = m*x + b
%             disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
%             PrfEcc2(depth_idx,hem) = c(1)*2 + c(2);
        end
    end
    %save([res_path 'pRFsizes.mat'],'PrfBinned','PrfEcc2','depth_val','Bins','hems')
    save([res_path 'pRFsizes.mat'],'Prfs','depth_val','Bins','hems')
    
end



