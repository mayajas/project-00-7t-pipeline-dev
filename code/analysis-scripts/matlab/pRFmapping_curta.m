function pRFmapping_curta(sub,depth_idx,ncpus)
if ncpus > 1
    p = parpool('local',ncpus);
    disp(['Parpool size: ' num2str(p.NumWorkers)])
end

project_name    = 'project-00-7t-pipeline-dev';

% set paths
scripts_path    = ['/home/mayaaj90/scripts/' project_name '/'];
data_path       = ['/scratch/mayaaj90/' project_name '/pRF/data/'];
FS_path         = ['/scratch/mayaaj90/' project_name '/pRF/data_FS/'];
programs_path   = '/home/mayaaj90/programs';
anat_path       = ['/scratch/mayaaj90/' project_name '/output/anat/_subject_id_'];

% add scripts path
addpath(scripts_path)

% add paths to required software
addpath(genpath(fullfile(programs_path,'vistasoft')))
addpath(genpath(fullfile(programs_path,'spm12')))
rmpath(genpath(fullfile(programs_path,'spm12','external','fieldtrip')))
addpath(genpath(fullfile(programs_path,'samsrf')))

% iterables
depth_text                       = {'-1.5','-1.0','-0.5','0.0','0.5','1.0','1.5','2.0','2.5'};

% get fMRI and pRF parameters
params = fMRIparameters(depth_text{depth_idx});
params = pRFparameters(params);                                      


% get the folder contents
thisSubject = ['sub-0' num2str(sub)];
surfsdir = [data_path 'surfs/_subject_id_' thisSubject];
maskdir = [data_path 'occLabels/_subject_id_' thisSubject];
fsdir       = [FS_path thisSubject];

% copy T1_out.nii to fsdir
if ~exist([fsdir filesep 'mri' filesep 'orig' filesep params.fMRI.struct.name '.nii'],'file')
    copyfile([anat_path thisSubject filesep params.fMRI.struct.name '.nii'],[fsdir filesep 'mri' filesep 'orig' filesep])
end

% run prf mapping
runPRFmapping(surfsdir,maskdir,fsdir,scripts_path,params,depth_text{depth_idx});

end                                            
                                            