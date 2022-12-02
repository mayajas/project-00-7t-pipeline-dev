function [params]  = pRFparameters(params)

% pRF parameters
params.pRF.fitHRF                = 0;        % whether to fit HRF using a run of full visual 
                                            % field stimulation data
                                            % if 0, then canonical HRF is used
if params.pRF.fitHRF
    params.pRF.HRF_name          = 'hrf';    % name of HRF run (ignored if fitHRF=0)
    params.pRF.HRF_stimdur       = 2;        % duriation of visual field stimulation block
                                            % (ignored if fitHRF=0)
end
params.pRF.num_stimuli           = 1;        % number of different types of pRF mapping stimuli used
                                            % here: 1 - classic(bar), 
                                            % 2 - task-relevant pRF mapping
params.pRF.stimulus_names        = {{'bar'}};
params.pRF.num_runs              = [{2}];    % number of runs per stimulus 
                                            % type
params.pRF.runs                  = [{1:2}];
params.pRF.max_ecc               = [{6.58}]; % maximum eccentricities per 
                                            % stimulus type
params.pRF.min_ecc               = [{0}];    % minimum eccentricities per 
                                            % stimulus type
params.pRF.only_positive         = {[0 0 1]}; % Which parameters must be positive? (x,y,sigma)
% params.pRF.ctxsteps              = {num2str(depth_val(depth_idx))}; % 'Cortex sampling steps?', 'Choosing default cortex sampling of 0.5.'
params.pRF.rulenum               = 7;   % RuleStrs = {'Mean' 'Median' 'Maximum' 'Minimum' 'Sum' 'Geomean' 'None'};
params.pRF.nrmls                 = 'Yes';    % 'Normalise time series?', '', 'Yes', 'No'
params.pRF.avrg                  = 'Average';% 'How to handle multiple runs?', '', 'Average', 'Concatenate', 'Separate'
params.pRF.TR                    = params.fMRI.EPI.TR;        % repetition time
end