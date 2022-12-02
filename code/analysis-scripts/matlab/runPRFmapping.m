function runPRFmapping(surfsdir,maskdir,fsdir,scripts_path,params,depth_text)    
    %%  Data organization for pRF model fitting
    if ~exist(fullfile([fsdir filesep 'fMRI'],[params.pRF.stimulus_names{1}{1}],['depth_' depth_text],['lh_surf_sess1.mgh'])) || ...
            ~exist(fullfile([fsdir filesep 'fMRI'],[params.pRF.stimulus_names{1}{1}],['depth_' depth_text],['rh_surf_sess1.mgh']))
        organizePreprocessedData(surfsdir, fsdir, params, depth_text);
    end
    

    %%   Create occipital labels
    if ~exist([fsdir filesep 'label' filesep 'lh_occ_depth' depth_text '.label'],'file') || ...
            ~exist([fsdir filesep 'label' filesep 'rh_occ_depth' depth_text '.label'],'file')
        if exist([maskdir filesep '_hemi_lh_sampling_range_' depth_text filesep 'lh_occ_depth' depth_text '.label'],'file') && ...
                exist([maskdir filesep '_hemi_rh_sampling_range_' depth_text filesep 'rh_occ_depth' depth_text '.label'],'file')
            
            copyfile([maskdir filesep '_hemi_lh_sampling_range_' depth_text filesep 'lh_occ_depth' depth_text '.label'],...
                [fsdir filesep 'label' filesep 'lh_occ_depth' depth_text '.label'])
            copyfile([maskdir filesep '_hemi_rh_sampling_range_' depth_text filesep 'rh_occ_depth' depth_text '.label'],...
                [fsdir filesep 'label' filesep 'rh_occ_depth' depth_text '.label'])
        end
    end

    %% Surface projection
    for pRF_stim = 1:params.pRF.num_stimuli
        for pRF_con = 1:length(params.pRF.stimulus_names{pRF_stim})
            stimulus_name = params.pRF.stimulus_names{pRF_stim}{pRF_con};
            
            %   Create Samsrf surface data files
            pRF_path = [fsdir filesep 'fMRI' filesep stimulus_name filesep 'depth_' depth_text filesep];
            params.pRF.path{pRF_stim}{pRF_con} = pRF_path;
            cd(params.pRF.path{pRF_stim}{pRF_con})
            
            if (~exist(['lh_surf_sess1.mat']) || ~exist(['rh_surf_sess1.mat']))
                hems = 1:2;
                for hem = hems
                    if hem == 1
                        hem_txt = 'lh';
                    elseif hem == 2
                        hem_txt = 'rh';
                    end

                    for run = 1:params.fMRI.EPI.num_runs(pRF_stim)
                        pRF_files{run,1} = [hem_txt '_surf_sess' num2str(run) '.mgh'];
                    end
                    SurfaceProjection(fsdir,pRF_files,params.pRF.path{pRF_stim}{pRF_con},params.pRF)
                end
            end
        end
    end

    %%   Preparing apertures
    for pRF_stim = 1:params.pRF.num_stimuli
        for pRF_con = 1:length(params.pRF.stimulus_names{pRF_stim})
            cd(params.pRF.path{pRF_stim}{pRF_con})
            if ~exist([params.pRF.path{pRF_stim}{pRF_con} 'aps_pRF.mat'])         
                delete([params.pRF.path{pRF_stim}{pRF_con} 'src_pRF_Gaussian.mat'])
                copyfile([scripts_path 'Apertures' filesep 'stimulus_' params.pRF.stimulus_names{pRF_stim}{pRF_con} '.mat'],[params.pRF.path{pRF_stim}{pRF_con} 'aps_pRF.mat'])
                ApFrm  = struct2cell(load('aps_pRF.mat'));
                ApFrm = ApFrm{1};
                save('aps_pRF.mat','ApFrm')
            end
        end
    end


    %%  RETINOTOPIC MAPPING
    for pRF_stim = 1:params.pRF.num_stimuli
        for pRF_con = 1:length(params.pRF.stimulus_names{pRF_stim})
            hems = 1:2;
            for hem = hems
                if hem == 1
                    hem_txt = 'lh';
                elseif hem == 2
                    hem_txt = 'rh';
                end

                %   Define the pRF model
                cd(params.pRF.path{pRF_stim}{pRF_con})
                if ~exist([hem_txt '_pRF_Gaussian.mat'],'file')
                    if exist([hem_txt '_surf.mat'])
                        SrfFiles = [hem_txt '_surf.mat'];
                    else
                        SrfFiles = [hem_txt '_surf_sess1.mat'];
                    end
                    
                    Roi = [fsdir filesep 'label' filesep hem_txt '_occ_depth' depth_text];

                    prf_params.TR = params.fMRI.EPI.TR;
                    prf_params.max_ecc = params.pRF.max_ecc{pRF_stim}(pRF_con);
                    prf_params.only_positive = params.pRF.only_positive{pRF_stim}(pRF_con,:);
                    prf_params.fitHRF = params.pRF.fitHRF;
                    if prf_params.fitHRF
                        prf_params.HRF_filepath = [params.fMRI.EPI.path{end} 'hrf_' hem_txt '_data4D.mat'];
                    end
                    
                    Standard_2d_Gaussian_Prf(params.pRF.path{pRF_stim}{pRF_con},SrfFiles, Roi, prf_params)
                end

%                 %   Post-processing pRF map
%                 if exist([hem_txt '_pRF_Gaussian.mat'],'file') && ~exist(['postproc_' hem_txt '_pRF_Gaussian.mat'],'file')
%                     load([hem_txt '_pRF_Gaussian.mat'],'Model','Srf')
%                     Roi = [fsdir filesep 'label' filesep hem_txt '_occ_depth' depth_text];
%                     Srf = samsrf_surfcalcs(Srf,Roi,0.01,[params.pRF.min_ecc{pRF_stim}(pRF_con) params.pRF.max_ecc{pRF_stim}(pRF_con)]);
%                     save(['postproc_' hem_txt '_pRF_Gaussian.mat'],'Model','Srf','-v7.3');
%                 end
                
%                 % Display maps
%                 %DisplayMaps
% 
% 
%                 % Delineate ROIs
%                 %DelineationTool
% 
%                 % Convert ROI labels to .nii
%                 if exist([params.pRF.path{pRF_stim}{pRF_con} 'ROIs_tproc_' hem_txt '_pRF_Gaussian'],'dir') 
%                     labels = dir(['ROIs_tproc_' hem_txt '_pRF_Gaussian' filesep '*.label']);
%                     for l = 1:length(labels)
%                         labelfile = [params.pRF.path{pRF_stim}{pRF_con} filesep 'ROIs_tproc_' hem_txt '_pRF_Gaussian' filesep labels(l).name];
%                         labelfile = labelfile(1:end-6); % remove extension
% 
%                         if ~exist([labelfile '.nii'],'file')
%                             funimg = ls([params.pRF.path{pRF_stim}{pRF_con} 'data4D*.nii']);
%                             funimg = [params.pRF.path{pRF_stim}{pRF_con} strtrim(funimg(1,1:end-4))];
%                             if ~isCross
%                                 strimg = [fsdir filesep  'mri' filesep  'orig' filesep  'longitudinalAnat'];
%                             else
%                                 strimg = [fsdir filesep  'mri' filesep  'orig' filesep  'anat'];
%                             end
%                             hemsurf = [fsdir filesep 'surf' filesep hem_txt];
%                             ctxsteps = 0.1;
%                             samsrf_label2nii(labelfile, funimg, strimg, hemsurf, ctxsteps, false, false)
%                         end
%                     end
%                 end

                       
           end         
        end
    end

end