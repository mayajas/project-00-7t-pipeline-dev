function [params] = fMRIparameters(depth)
% general
params.fMRI.B0                   = 7;        % magnetic field strength (Tesla)

% fMRI
params.fMRI.EPI.num_runs         = 2;        % number of runs per stimulus type
params.fMRI.EPI.TR               = 3;        % repetition time

bar1_subdir = '_sess_id_task-bar_run-01_sess_nr_0_sess_nvol_124';
bar2_subdir = '_sess_id_task-bar_run-02_sess_nr_1_sess_nvol_124';
lh_subdir   = '_hemi_lh_sampling_range_';
rh_subdir   = '_hemi_rh_sampling_range_';

params.fMRI.EPI.fMRI_name        = {[bar1_subdir filesep lh_subdir depth filesep 'lh.atask-bar_run-01_roi_warp4D_trans.mgh'],...
                                    [bar2_subdir filesep lh_subdir depth filesep 'lh.atask-bar_run-02_roi_warp4D_trans.mgh'],... 
                                    [bar1_subdir filesep rh_subdir depth filesep 'rh.atask-bar_run-01_roi_warp4D_trans.mgh'],...
                                    [bar2_subdir filesep rh_subdir depth filesep 'rh.atask-bar_run-02_roi_warp4D_trans.mgh']};

% params.fMRI.EPI.fMRI_name        = {['lh_bar_sess1_depth' depth '.mgh'],['lh_bar_sess2_depth' depth '.mgh'],...
%                                         ['rh_bar_sess1_depth' depth '.mgh'],['rh_bar_sess2_depth' depth '.mgh']};
params.fMRI.EPI.total_num_runs   = sum(params.fMRI.EPI.num_runs);

% struct
params.fMRI.struct.name          = 'UNI_corrected'; % anatomical image name
params.fMRI.struct.isMT          = 0;        % is the anatomical image an MT map
                                            % if not, it's assumed to be a T1 mprage
end