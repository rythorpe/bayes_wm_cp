clear

% Goals:
% 1) Compare simulated trial outcomes between the Reduced Bayesian and
%    WM Models at different WM buffer sizes.
% 2) Plot example outcomes from the two different models and the distance
%    between them.
%
% Implementational note: this script uses the old data format where all
% trials across subjects are concatenated along a single dimension.

% add paths
whichComp=3;
switch whichComp
    case 1 % jinrui
        basePath='/Users/apple/Dropbox/';
    case 2 % matt
        basePath= '~/Dropbox';
    case 3 % ryan
        basePath= '~/bayes_wm_cp';
    case 4 % mara
        basePath = '/Users/elenaoancea/Dropbox';
end

% addpath(genpath(fullfile(basePath, 'sharedMatlabUtilities')));
% addpath(genpath(fullfile(basePath, 'onlineCannonTask/modelFitting')));
addpath(genpath(fullfile(basePath)));

%% Create synthetic outcome data %%

numBlocks = 3;
trials_per_block = 1000;
num_trials = trials_per_block * numBlocks;
colors = 1:3;
blocks = 1:numBlocks;
setSizes = [1,2,3]; % corresponds to different blocks
haz = 0.1;
noiseStd = deg2rad(20);
trial_outcome = nan(num_trials, 1);
trial_mean = nan(num_trials, 1);
trial_cp = zeros(num_trials, 1);
trial_runlength = zeros(num_trials, 1);
trial_set_size = nan(num_trials, 1);
trial_color = nan(num_trials, 1);
trial_block = nan(num_trials, 1);
blkColBeg = zeros(num_trials, 1);

rng('default')
for block_idx=1:numBlocks
    numColors = block_idx;
    trials_per_color = trials_per_block/numColors;
    for color_idx=1:numColors
        % first trial for the first color occurs at the start of the
        % block, and will reoccur after the other colors have been
        % sampled round-robin (1,2,3,1,2,3...)
        start_idx = (block_idx-1)*trials_per_block + ...
                    color_idx;
        end_idx = (block_idx-1)*trials_per_block + ...
                   trials_per_block;
        % note beginning index of color block
        blkColBeg(start_idx) = 1;
        % set initial mean and iterate over trials for this color
        current_mean = rand*2*pi;
        runlength = 0;
        for trial_idx=start_idx:numColors:end_idx
            if rand <= haz
                % update mean at random
                current_mean = rand*2*pi;
                trial_cp(trial_idx) = 1;
                runlength = 0;
            end
            trial_runlength(trial_idx) = runlength;
            trial_mean(trial_idx) = current_mean;
            trial_outcome(trial_idx) = circ_vmrnd(current_mean, 1/noiseStd^2, 1);
            trial_set_size(trial_idx) = numColors;
            trial_color(trial_idx) = colors(color_idx);
            trial_block(trial_idx) = blocks(block_idx);
            runlength = runlength + 1;
        end
    end
end

%% Simulate predicted outcomes for synthetic subjects: WM model %%

% [buffSize,haz,motorStd,cpLikeBias,noiseStd]
% input.startPoint = [1, 0.1, pi/8, 0.5, deg2rad(20)]; % Note: doesn't matter much with data.searchStartPoints=true
input.ub = [8, 0.5, pi, 1, pi];
input.lb = [0, 0.001, pi/180, 0.01, pi/180];
input.whichParams = [1, 1, 1, 1, 1];
input.coordSys = 'polar';

% Set ground-truth values for parameters
numParams = length(input.ub);
buff_sizes = 1:0.2:15;
motorNoiseStd = 1e-3; % careful!!! don't make too small or circ_vmrnd will get stuck in infinite while-loop
cpLikeBias = 1;

numBuffSizes = length(buff_sizes);
mod_params_wm = nan(numBuffSizes,numParams);
for bi=1:numBuffSizes
    mod_params_wm(bi,:) = [buff_sizes(bi), haz, motorNoiseStd, cpLikeBias, noiseStd];
end

% Fitting/simulation settings
input.doFit = false;
input.searchStartPoints = false;
input.doSim = true;
input.usePrior = false;

% Generate data
trial_prediction_wm = nan(length(trial_outcome),numBuffSizes);
for bi=1:numBuffSizes
    % set params to ground-truth values for simulation
    input.startPoint = mod_params_wm(bi,:);
    
    % select input data
    input.beliefs = []; % Note: doesn't matter with input.doFit=false
    input.outcomes = trial_outcome;
    input.blocks = trial_block;
    input.setSizes = trial_set_size;
    input.colors = trial_color;
    
    % simulate model predictions
    output = fit_AL_WM_allSetSize(input);

    % store simulated belief data
    trial_prediction_wm(:,bi) = output.simBelief;
end
clear input

%% Simulate predicted outcomes for synthetic subjects: reduced Bayesian model %%
% 1) haz
% 2) likelihood weight
% 3) update standard deviation intercept
% 4) update standard deviation slope (over absolute error magnitude)
% 5) max prob of binary updates
% 6) width of binary update distribution
% 7) memory capacity
update_std_slope = 0; % no adaptive "motor noise" based on error-driven update
max_bin_update = 0; % set chance of no update to 0
mem_capacity = 4;
mod_params_rb = [haz, cpLikeBias, motorNoiseStd, 0, max_bin_update, 4, mem_capacity];

input = struct;
input.lb = [10e-10, 10e-10, .0001, -1, 0, 0, 0];
input.ub = [1-10e-10, 1-10e-10, 4, 1, 1, 4, 3];
input.startPoint = mod_params_rb;
input.whichParams = [true(1,7)];
input.doFit = false;
input.usePrior = false;
input.simFromScratch = true;

% note: fitRedBayesModel_zombies_allSetSize.m has now been modified to
% simulate data for multiple set sizes (each in a different block) all at
% once
input.PE = ones(size(trial_outcome)) .* 0.5; % arbitary for simulated model data
input.UP = ones(size(trial_outcome)); % arbitary for simulated model data
input.skipLike = true(size(trial_outcome));
input.nn = ones(size(trial_outcome)).*noiseStd;
input.newBlock = blkColBeg;
input.outcomes = trial_outcome;
input.color = trial_color;
input.setSize = trial_set_size;

% simulate model predictions
output = fitRedBayesModel_zombies_allSetSize(input);            
trial_prediction_rb = output.modPred;

%% plot example trials from task %%
% compare reduced Bayesian predictions to WM model (with buff_size=1, i.e., LR of 1) predictions
sel_buffSize = buff_sizes==1;
plot_task_multicolor(trial_color, ...
                     trial_mean, ...
                     trial_outcome, ...
                     trial_prediction_rb, ... % subj's belief = redBayes model prediction
                     false(size(trial_outcome)), ... % no outlier trials; all are false
                     trial_prediction_wm(:,sel_buffSize), ... % WM model prediction
                     trial_cp, ...
                     [300, 400], ... % select bounds of a group of trials for plotting
                     'polar3D')

%% Compare distance between model predictions %%
figure
legend_labels = cell(size(setSizes));
hold on
for ssi=1:length(setSizes)
    sel = trial_set_size==setSizes(ssi);
    model_diffs = circ_dist(repmat(trial_prediction_rb(sel),1,numBuffSizes), trial_prediction_wm(sel,:));
    plot(buff_sizes, mean(abs(model_diffs),1))
    legend_labels{ssi} = sprintf('set size %d',setSizes(ssi));
end
hold off
legend(legend_labels)
xlabel('buffer size')
ylabel('mean distance between model predictions')
