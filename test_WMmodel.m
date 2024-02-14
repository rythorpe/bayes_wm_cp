clear

% Goals:
% 1) Conduct a parameter recovery test for the working memory model
%    (fit_AL_WM_allSetSize.m).
% 2) Demo trial-after-changepoint (TAC) analysis.
% 3) Plot learning-rate (LR) curves for different levels of likelihood
%    weight to demo how it reduces the LR slope (also depends on polar vs.
%    non-corrected polar LR).

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

numSubjs = 30; % try 8 for shorter runtime
numBlocks = 3;
trials_per_block = 300;
trials_per_subj = numBlocks*trials_per_block;
colors = 1:3;
blocks = 1:3;
haz = 0.2;
noiseStd = deg2rad(20);
subjList = 1:numSubjs;
data = struct([]);

rng('default')
for subj_idx=1:numSubjs
    data(subj_idx).prediction = nan(trials_per_subj, 1); % for synthetic belief data: remains empty until later
    data(subj_idx).outcome = nan(trials_per_subj, 1);
    data(subj_idx).mean = nan(trials_per_subj, 1);
    data(subj_idx).cp = zeros(trials_per_subj, 1);
    data(subj_idx).runlength = zeros(trials_per_subj, 1);
    data(subj_idx).set_size = nan(trials_per_subj, 1);
    data(subj_idx).color = nan(trials_per_subj, 1);
    data(subj_idx).block = nan(trials_per_subj, 1);
    data(subj_idx).subj = subj_idx + zeros(trials_per_subj, 1); % only used when concatenating trials across subjects near end of script
    for block_idx=1:numBlocks
        numColors = block_idx;
        trials_per_color = trials_per_block/numColors;
        for color_idx=1:numColors
            % first trial for the first color occurs at the start of the
            % block, and will reoccur after the other colors have been
            % sampled round-robin (1,2,3,1,2,3...)
            start_idx = (block_idx-1)*trials_per_block + color_idx;
            end_idx = (block_idx-1)*trials_per_block + trials_per_block;
            % set initial mean and iterate over trials for this color
            current_mean = rand*2*pi;
            runlength = 0;
            for trial_idx=start_idx:numColors:end_idx
                if rand <= haz
                    % update mean at random
                    current_mean = rand*2*pi;
                    data(subj_idx).cp(trial_idx) = 1;
                    runlength = 0;
                end
                data(subj_idx).runlength(trial_idx) = runlength;
                data(subj_idx).mean(trial_idx) = current_mean;
                data(subj_idx).outcome(trial_idx) = circ_vmrnd(current_mean, 1/noiseStd^2, 1);
                data(subj_idx).set_size(trial_idx) = numColors;
                data(subj_idx).color(trial_idx) = colors(color_idx);
                data(subj_idx).block(trial_idx) = blocks(block_idx);
                runlength = runlength + 1;
            end
        end
    end
end

%% Simulate predicted outcomes for synthetic subjects%%

% [buffSize,haz,motorStd,cpLikeBias,noiseStd]
% input.startPoint = [1, 0.1, pi/8, 0.5, deg2rad(20)];
input.ub = [8, 0.5, pi, 1, pi];
input.lb = [0, 0.001, pi/180, 0.01, pi/180];
input.whichParams = [1, 1, 1, 1, 1];
input.coordSys = 'polar';

% set ground-truth values for parameters
numParams = length(input.ub);
pop_param_mean = [3.0, 0.3, deg2rad(10), 0.7, deg2rad(15)];
pop_param_std = [0.8, 0.1, deg2rad(3), 0.2, deg2rad(4)];
params_true = nan(numSubjs,numParams);
for si=1:numSubjs
    sampled_params = normrnd(pop_param_mean, pop_param_std);
    % make sure none of the params go outside bounds
    for param_idx=1:numel(sampled_params)
        ub = input.ub(param_idx);
        lb = input.lb(param_idx);
        % if so, resample
        while sampled_params(param_idx) < lb || sampled_params(param_idx) > ub
            sampled_params(param_idx) = normrnd(pop_param_mean(param_idx), pop_param_std(param_idx));
        end
    end
    params_true(si,:) = sampled_params;
end
% pop_ub = [6.0, 0.50, deg2rad(15), 1, pi/4];
% pop_lb = [0.0, 0.001, deg2rad(1), 0.01, deg2rad(1)];
% params_true = rand(numSubjs,numParams) .* repmat(pop_ub-pop_lb,numSubjs,1) + repmat(pop_lb,numSubjs,1);

% fitting/simulation settings
input.doFit = false;
input.searchStartPoints = false;
input.doSim = true;
input.usePrior = false;

% generate data
for si=1:numSubjs
    % set params to ground-truth values for simulation
    input.startPoint = params_true(si,:);
    
    % select input data
    input.beliefs = []; % Note: doesn't matter with input.doFit=false
    input.outcomes = data(si).outcome;
    input.blocks = data(si).block;
    input.setSizes = data(si).set_size;
    input.colors = data(si).color;
    
    % run simulation
    output = fit_AL_WM_allSetSize(input);

    % store simulated belief data
    data(si).prediction = output.simBelief;
end
clear input

%% Plot example trials from task

% select a subject; belief from subjects with parameters closer to an ideal observer will be easier to interpret
si = 2;
plot_task_multicolor(data(si).color, ...
                     data(si).mean, ...
                     data(si).outcome, ...
                     data(si).prediction, ... % for now we'll assume the observer's belief matches the model's belief
                     zeros(trials_per_subj,1), ... % no outlier trials; all are false
                     data(si).prediction, ... % model's belief
                     data(si).cp, ...
                     [625, 675], ... % select bounds of a group of trials for plotting
                     'polar3D')

%% Test parameter recovery by then fitting synthetic data %%

% fitting/simulation settings
input.doFit = true;
input.searchStartPoints = true;
input.doSim = false;
input.usePrior = false;

% [buffSize,haz,motorStd,cpLikeBias,noiseStd]
input.startPoint = [1, 0.1, deg2rad(10), 0.7, deg2rad(20)];
% fit with the following bounds (note that this was already set above)
input.ub = [8, 0.5, pi, 1, pi];
input.lb = [0, 0.001, pi/180, 0.01, pi/180];
input.whichParams = [1, 1, 1, 1, 1];
input.coordSys = 'polar';

% do first round of parameter fitting to estimate prior
params_fit_recov = nan(numSubjs,numParams);
fit_times = nan(1,numSubjs);
for si=1:numSubjs

    % select input data
    input.beliefs = data(si).prediction;
    input.outcomes = data(si).outcome;
    input.blocks = data(si).block;
    input.setSizes = data(si).set_size;
    input.colors = data(si).color;
    
    % run fit and log elapsed time for each iteration
    tic
    output = fit_AL_WM_allSetSize(input);
    fit_times(si) = toc;
    params_fit_recov(si,:) = output.params; % Store optimized params
    disp(['Subj. #',num2str(si),': ',num2str(output.params)])
end

mean(fit_times)  % 41.0s per subj w/parfor; 148.6s w/for

% expectation maximization (with estimated prior)
numIt = 2;
for i=1:numIt
    
    % now run optimization using the prior
    input.doFit = true;
    input.searchStartPoints = true;
    input.doSim = false;
    input.usePrior = true;
    
    % define group-level stats
    input.paramMean = mean(params_fit_recov);
    input.paramCov = cov(params_fit_recov);
    
    % iterate through subjects
    for si=1:numSubjs
        
        % initialize to the previously-determined param values
        % ...doesn't matter much since we are searching start points
        input.startPoint = params_fit_recov(si,:);
        
        % select input data
        input.beliefs = data(si).prediction;
        input.outcomes = data(si).outcome;
        input.blocks = data(si).block;
        input.setSizes = data(si).set_size;
        input.colors = data(si).color;
        
        % run fit
        output = fit_AL_WM_allSetSize(input);
        params_fit_recov(si,:) = output.params; % Store optimized params
        disp(['Subj. #',num2str(si),': ',num2str(output.params)])
    end
end

%% Plot results of parameter recovery experiment %%

param_names = {'buffer size','hazard rate','motor std','CP bias','noise std'};

% plot ground-truth param histograms
figure
for pari=1:numParams
    param_mean = mean(params_true(:,pari));
    subplot(1,numParams,pari)
    histogram(params_true(:,pari),15)
    hold on
    line([param_mean param_mean],ylim,'color','r')
    hold off
    xlabel(param_names{pari})
    ylabel('counts')
end

% plot recovered param histograms
figure
for pari=1:numParams
    param_mean = mean(params_fit_recov(:,pari));
    subplot(1,numParams,pari)
    histogram(params_fit_recov(:,pari),15)
    hold on
    line([param_mean param_mean],ylim,'color','r')
    hold off
    xlabel(param_names{pari})
    ylabel('counts')
end

% plot parameter recovery correlations
predictor = params_true;
response = params_fit_recov;
plot_paramFitRecovery(predictor,response,param_names)

%% Plot trials-after-changepoint analysis %%

LR_polar = num2cell(nan(numSubjs, 3));
LR_noCorrect = num2cell(nan(numSubjs, 3));

for subj_i=1:numSubjs
    for ss_i=1:3
        subj_data = data(subj_i);
        setSize = blocks(ss_i); % note: set size happens to be the same as the block label, but this doesn't necessarily need to be true
        sel = subj_data.set_size==setSize;
        allCols = unique(subj_data.color(sel));
        blkColBeg = [1; zeros(trials_per_block/setSize, 1)];
        TAC = subj_data.runlength(sel) + 1; % double-check: runlength of 0 == TAC of 1
        
        UP_noCorrect = nan(sum(sel), 1);
        PE_noCorrect = nan(sum(sel), 1);
        UP_polar = nan(sum(sel), 1);
        PE_polar = nan(sum(sel), 1);
        for col_idx=1:length(allCols) % loop through colors
            colSel = sel & subj_data.color==allCols(col_idx);
            colSel_in_ss = subj_data.color(sel)==allCols(col_idx);
            [~,UP_noCorrect(colSel_in_ss),PE_noCorrect(colSel_in_ss)] = computeLR(subj_data.outcome(colSel),subj_data.prediction(colSel),blkColBeg,'polarNoCorrect');
            [~,UP_polar(colSel_in_ss),PE_polar(colSel_in_ss)] = computeLR(subj_data.outcome(colSel),subj_data.prediction(colSel),blkColBeg,'polar');
        end
        
        xMat = [ones(sum(sel), 1), PE_noCorrect.*(TAC==1), PE_noCorrect.*(TAC==2), PE_noCorrect.*(TAC==3), PE_noCorrect.*(TAC>=4)];
        yVar = UP_noCorrect;
        LR_noCorrect{subj_i,ss_i,:} = regress(yVar, xMat);
        LR_noCorrect{subj_i,ss_i} = reshape(LR_noCorrect{subj_i,ss_i}', [], 5);
        
        xMat = [ones(sum(sel),1), PE_polar.*(TAC==1), PE_polar.*(TAC==2), PE_polar.*(TAC==3), PE_polar.*(TAC>=4)];
        yVar = UP_polar;
        LR_polar{subj_i,ss_i,:} = regress(yVar, xMat);
        LR_polar{subj_i,ss_i} = reshape(LR_polar{subj_i,ss_i}', [], 5);
    end
end

TACanalysis(LR_noCorrect,LR_polar,true(numSubjs,1));

%% Plot learning-rate curves for different levels of likelihood weight (cpLikeBias) %%

likeWeights = [1.0, 0.75, 0.5, 0.25, 0.0]; % number of different likelihood weights to try

% [buffSize,haz,motorStd,cpLikeBias,noiseStd]
% input.startPoint = [1, 0.1, pi/8, 0.5, deg2rad(20)];
input.ub = [8, 0.5, pi, 1, pi]; % Note: doesn't matter much with input.doFit=false
input.lb = [0, 0.001, pi/180, 0.01, pi/180]; % Note: doesn't matter much with doFit=false
input.whichParams = [1, 1, 1, 1, 1];
input.coordSys = 'polar';

% fitting/simulation settings
input.doFit = false;
input.searchStartPoints = false;
input.doSim = true;
input.usePrior = false;

% generate data
outcomes_lw = nan(numSubjs*trials_per_block,length(likeWeights));
predictions_lw = nan(numSubjs*trials_per_block,length(likeWeights));
blkBeg = repmat([1; zeros(trials_per_block, 1)],[numSubjs, length(likeWeights)]);
legend_labels = cell(size(likeWeights));
for wi=1:length(likeWeights)
    % set params to for simulation
    % use small motor noise to make learning rate trend clearer
    input.startPoint = [5.0, haz, 0.001, likeWeights(wi), noiseStd];
    
    % select input data
    sel = cat(1, data.set_size)==1; % use set size 1 trials across all subjects
    all_outcomes = cat(1, data.outcome);
    all_blocks = cat(1, data.subj); % in this case blocks belong to different subjects
    all_setSizes = cat(1, data.set_size);
    all_colors = cat(1, data.color);

    input.beliefs = []; % Note: doesn't matter with input.doFit=false
    input.outcomes = all_outcomes(sel);
    input.blocks = all_blocks(sel);
    input.setSizes = all_setSizes(sel);
    input.colors = all_colors(sel);
    
    % run simulation
    output = fit_AL_WM_allSetSize(input);

    % store simulated belief data
    predictions_lw(:,wi) = output.simBelief;
    outcomes_lw(:,wi) = all_outcomes(sel);

    % set legend name for current likelihood weight
    legend_labels{wi} = sprintf('%0.2f',likeWeights(wi));
end

% plot LR as a function of prediction error
% note: using the average of the polar and non-corrected polar LR (default)
% has an interesting effect that helps us visualize what's going on
% see plot_LR.m function for more details
plot_LR(outcomes_lw,predictions_lw,blkBeg)
legend(legend_labels)
