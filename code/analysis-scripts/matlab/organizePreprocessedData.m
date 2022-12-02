function [databasedir] = organizePreprocessedData(databasedir, fsdir, params, depth_text)

disp('Now organizing the preprocessed data...')
savedir = [fsdir filesep 'fMRI'];
if ~exist(fullfile(savedir))
    mkdir(fullfile(savedir))
end

% hrf
if params.pRF.fitHRF
    hrf = dir(fullfile(databasedir,params.fMRI.EPI.fMRI_name{3}));
    hrf = hrf(3:end);
    for i = 1:length(hrf)
        files = dir(fullfile(databasedir,params.fMRI.EPI.fMRI_name{3},hrf(i).name,'rbuafPR*'));
        nii = niftiRead(fullfile(databasedir,params.fMRI.EPI.fMRI_name{3},hrf(i).name,files(11).name));
        data4D = nii.data;
        for k=11:length(files)-1 - 2
            nii = niftiRead(fullfile(databasedir,params.fMRI.EPI.fMRI_name{3},hrf(i).name,files(k+1).name));
            data4D = cat(4, data4D, nii.data);
        end
        if ~exist(fullfile(savedir,'hrf'))
            mkdir(fullfile(savedir,'hrf'))
        end
        dtiWriteNiftiWrapper(data4D, nii.qto_xyz,fullfile(savedir,'hrf','data4D.nii') , 1, '', [],[],[],[], 1);
        gunzip(fullfile(savedir,'hrf','*.gz'))
    end
end

% pRF mapping functional runs
for stim = 1:length(params.fMRI.EPI.fMRI_name)
    pRF_sess = dir(fullfile(databasedir,params.fMRI.EPI.fMRI_name{stim}));
    pRF_file   = fullfile(databasedir,params.fMRI.EPI.fMRI_name{stim});
    
    if length(pRF_sess) ~= 1%sum(params.pRF.num_runs{stim})*2 %lh and rh
        error('Number of pRF sessions doesn''t match design')
    end

    % get session number, hemisphere
    sess_nr = num2str(str2double(extractBefore(extractAfter(pRF_sess.name,'run-'),3)));
    if contains(pRF_sess.name,'lh') && ~contains(pRF_sess.name,'rh')
        hem = 'lh';
    elseif contains(pRF_sess.name,'rh') && ~contains(pRF_sess.name,'lh')
        hem = 'rh';
    else
        error('Can''t tell if this is right or left hemisphere!')
    end
    
    % verify that depth matches the intended depth
    if ~contains(pRF_sess.folder,depth_text) 
        error('Can''t tell if this is the right surface projection depth')
    end
    
    % verify that the stimulus type matches the intended one
    stimulus_name = params.pRF.stimulus_names{1}{1};
    if ~contains(pRF_sess.folder,stimulus_name) || ~contains(pRF_sess.name,stimulus_name)
        error(['Check that this file corresponds to the ' stimulus_name ' stimulus.'])
    end

    % make dir if doesn't yet exist
    if ~exist(fullfile(savedir,[stimulus_name],['depth_' depth_text]))
        mkdir(fullfile(savedir,[stimulus_name],['depth_' depth_text]))
    end

    % copy functional files 
    if ~exist(fullfile(savedir,[stimulus_name],['depth_' depth_text],[hem '_surf_sess' sess_nr '.mgh']))
        copyfile(pRF_file,fullfile(savedir,[stimulus_name],['depth_' depth_text],[hem '_surf_sess' sess_nr '.mgh']))
    else
        display(['File ' [hem '_surf_sess' sess_nr '.mgh'] ' already exists for depth ' depth_text])
    end
    
end

