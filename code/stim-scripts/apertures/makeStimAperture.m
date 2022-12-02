close all
clear all;

behav_dir = '/home/mayajas/Documents/project-00-7t-pipeline-dev/data/raw/sub-01/behav/';
stim_dir = '/home/mayajas/Documents/project-00-7t-pipeline-dev/code/stim-scripts/apertures/';
mkdir([stim_dir 'bar_run-01'])
mkdir([stim_dir 'bar_run-02'])


% load([behav_dir 'data20191021T105933.mat'])
% load([behav_dir 'data20191021T112726.mat'])
load([behav_dir 'data20191021T114252.mat'])
orig_size = size(original_stimulus.images{1});
img = original_stimulus.images{1};
stim_size = 200;
% stim = zeros(orig_size(1),orig_size(2),...
%     size(stimulus.seq,1)/60);
stim = zeros(stim_size,stim_size,...
    size(stimulus.seq,1)/60);
size(stim)
%%
figure
i = 1;
for f = 1:60:size(stimulus.seq,1)
  c = squeeze(img(:,:,stimulus.seq(f)));
  idx = zeros(size(c));
  idx(c==1 | c==254) = 1;
  stim(:,:,i) = imresize(idx, [200 200], 'nearest');
  imagesc(stim(:,:,i)); colorbar;
  pause(0.01);
  %f
  i = i + 1;
end

% stim = stim(:,:,1:60:size(stim,3));
stim = stim(:,:,5:end);
figure
for f = 1:size(stim,3)
  imagesc(stim(:,:,f)); colorbar;
  pause(0.5);
  imwrite(stim(:,:,f),[stim_dir 'bar_run-01' filesep 'frame_00' num2str(f) '.png'])
  imwrite(stim(:,:,f),[stim_dir 'bar_run-02' filesep 'frame_00' num2str(f) '.png'])
end

save('stimulus_bar','stim');