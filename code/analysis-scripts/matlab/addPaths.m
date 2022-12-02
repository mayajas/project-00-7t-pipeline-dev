function [scripts_path, data_path] = addPaths()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add all paths needed for CVPL_main.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = warning ('off','all');

main_path=fileparts(fileparts(matlab.desktop.editor.getActiveFilename));

% add scripts path
addpath(genpath(fullfile(main_path,'scripts')))

% script, data and data_FS path
scripts_path = [fullfile(main_path,'scripts') filesep];
data_path = [fullfile(main_path,'data') filesep];

% add paths to required software
programs_path = '/home/mayajas/Documents/programs';
addpath(genpath(fullfile(programs_path,'vistasoft')))
addpath(genpath(fullfile(programs_path,'spm12')))
rmpath(genpath(fullfile(programs_path,'spm12','external','fieldtrip')))
addpath(genpath(fullfile(programs_path,'samsrf')))
addpath(genpath(fullfile(programs_path,'LpsyToolbox')))
% addpath(genpath(fullfile(programs_path,'PsychToolbox')))


end