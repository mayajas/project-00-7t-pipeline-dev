function pRFmapping_main()
dbstop if error

project_name    = 'project-00-7t-pipeline-dev';

scripts_path    = ['/home/mayajas/Documents/' project_name '/code/analysis-scripts/matlab/'];
data_path       = ['/home/mayajas/scratch/' project_name '/pRF/data/'];
FS_path         = ['/home/mayajas/scratch/' project_name '/pRF/data_FS/'];
programs_path   = '/home/mayajas/Documents/programs';
anat_path       = ['/home/mayajas/scratch/' project_name '/output/anat/_subject_id_'];

% add scripts path
addpath(scripts_path)

% add paths to required software
addpath(genpath(fullfile(programs_path,'vistasoft')))
addpath(genpath(fullfile(programs_path,'spm12')))
rmpath(genpath(fullfile(programs_path,'spm12','external','fieldtrip')))
addpath(genpath(fullfile(programs_path,'samsrf_local')))

% iterables
depth_text                       = {'-1.5','-1.0','-0.5','0.0','0.5','1.0','1.5','2.0','2.5'};

for sub = 2%1:4
    thisSubject = ['sub-0' num2str(sub)];
    surfsdir = [data_path 'surfs/_subject_id_' thisSubject];
    maskdir = [data_path 'occLabels/_subject_id_' thisSubject];
    fsdir       = [FS_path thisSubject];
    
    for depth_idx = 4%1:length(depth_text)
        % get fMRI and pRF parameters  
        params = fMRIparameters(depth_text{depth_idx});
        params = pRFparameters(params); 
        
        % copy T1_out.nii to fsdir
        if ~exist([fsdir filesep 'mri' filesep 'orig' filesep params.fMRI.struct.name '.nii'],'file')
            copyfile([anat_path thisSubject filesep params.fMRI.struct.name '.nii'],[fsdir filesep 'mri' filesep 'orig' filesep])
        end
        
        runPRFmapping(surfsdir,maskdir,fsdir,scripts_path,params,depth_text{depth_idx});
    end
end
          
end
                                            