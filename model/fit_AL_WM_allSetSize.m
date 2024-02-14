function [output]=fit_AL_WM_allSetSize(input)
% GOAL:
% Fit the adaptive learning working memory model to change point data with
% multiple set sizes
%
% Strategy:
% 1) make input include a vector of "set sizes" or a list of set sizes for
%    each block. 
% 2) loop through blocks when computing likelihood
% 3) divide effective memory capacity by set size in each block
% 4) loop within colors per block
% 5) sum log likelihood across colors and blocks.
%
% Take structure of inputs and turn it into a structure out outputs
%
% Inputs:
% -------
% .outcomes          -- what outcomes did subject see?
% .beliefs           -- where did participant put bucket on each trial?
% .blocks            -- vector of n categorical values denoting block 0 or block 
%                       0 for n trials
% .setSizes          -- vector of set sizes associated with each trial
%                       (should be constant across a given block)
% .colors            -- vector of categorical color values associated with
%                       each trial
% .ub                -- upper bound for each parameter
% .lb                -- lower bound for each parameter
% .startPoint        -- starting point for the search 
% .whichParams       -- vector of logicals saying which things will be fit.
%                       1) buffer size
%                       2) hazard
%                       3) motor standard deviation
%                       4) CP likelihood bias (exponent in range [0,1])
%                       5) noise block 1
%                       6) noise block 2
% .paramMean         -- 1-by-6 matrix of group-mean parameter values to use
%                       when computing the prior (only relevant when 
%                       .usePrior=true)
% .paramCov          -- 6-by-6 covariance matrix of group-level parameter 
%                       values to use when computing the prior (only 
%                       relevant when .usePrior=true)
% .doSim             -- logical value; if true, we'll simulate some data 
%                       from the best fitting model 
% .doFit             -- logical value; if true, fit model parameters. If 
%                       false, just compute local likelihood.
% .searchStartPoints -- logical value; if true, search among a uniform
%                       distribution of start points for the best fitting 
%                       model parameters. If false, use .startPoint
% .usePrior          -- logical value; if true, incorporate a normal prior
%                       distribution into the calculation of the posterior
% .coordSys          -- string identifying the coordinate system in which
%                       the task is performed (i.e., 'cartesian' or
%                       'polar'). Defaults to 'cartesian'
%
% Outputs:
% --------
% .params            -- final parameters
% .logLike           -- log likelihood of belief data given the input 
%                       parameters or best fitting model parameters
% .logLikeVec        -- log likelihood for each trial
% .simBelief         -- beliefs that are simulated from the model with best
%                       fitting parameters

% User parameter!!!
% Number of additional, randomly selected start points to try if
% optimization involves searching among various initial conditions
numSP = 15;

% Determine block boundaries
begBlock = [1; diff(input.blocks)~=0];
begBlockInd = find(begBlock==1);
endBlockInd = [begBlockInd(2:end)-1; length(input.blocks)];
blockTypes = unique(input.blocks);

% Vector of logicals indicating which parameters to optimize over during
whichParams = logical(input.whichParams);

% Initialize optional parameters if they haven't already been defined in
% the 'input' struct
if ~isfield(input, 'whichParams')
    input.whichParams = [1, 1, 1, 1, 1, 1];
end

if ~isfield(input, 'doFit')
    input.doFit = true;
end

if ~isfield(input, 'doSim')
    input.doSim = false;
end

if ~isfield(input, 'searchStartPoints')
    input.searchStartPoints = false;
end

if ~isfield(input, 'usePrior')
    input.usePrior = false;
end

if ~isfield(input, 'coordSys')
    input.coordSys = 'cartesian';
end

% Check for NaNs in the observer's beliefs (predictions)
if sum(isnan(input.beliefs)) + sum(isnan(input.outcomes)) > 0
    error('at least one trial has an invalid belief or outcome of NaN')
end

% Check that set size should matches number of unique colors for each block
for b_i=1:numel(endBlockInd)
    blockType_ind = find(blockTypes==input.blocks(begBlockInd(b_i)));
    blSel=input.blocks==blockTypes(blockType_ind);
    blockSetSize=unique(input.setSizes(blSel)); % Divide buffer size by number of colors!!!
    blockColors=unique(input.colors(blSel)); % Divide buffer size by number of colors!!!
    if length(blockSetSize) ~= 1
        error(['there should only be one set size per block, got ' ...
               '%d different sizes'], length(blockSetSize))
    end
    if blockSetSize ~= length(blockColors)
        error(['the set size for this block should reflect the ' ...
               'number of unique colors, got set size %d and %d ' ...
               'unique colors'], blockSetSize, length(blockColors))
    end
end

% Check that lower and upper bounds inputs are congruent
numParams = length(whichParams);
if numParams ~= length(input.lb) || numParams ~= length(input.ub)
    error(['lower bounds and upper bounds must each be a vector with ' ...
           'the same length as whichParams (%d)'], numParams)
end

% Is the noise param unique for each block type?
if numParams - 4 ~= length(blockTypes)
    uniqueNoiseBlocks = false;
else
    uniqueNoiseBlocks = true;
end

if input.searchStartPoints==false
    numSP = 1;
end
% initialize start point (universal for fitting, simulating, or just
% computing NLL)
startPoint = input.startPoint;

% Run optimization to fit parameters
% options = optimset('Algorithm','interior-point', 'MaxPCGIter', 5000, 'MaxProjCGIter', 5000, 'MaxSQPIter', 5000);
options = optimoptions(@fmincon,'display','notify');
numParams_opt = sum(whichParams);
if input.doFit

    % Validate congruency of input data sizes
    if length(input.outcomes)~=length(input.beliefs)
        warning("incongruent number of beliefs and outcomes!!")
    end

    % Set bounds for -log(likelihood) minimization
    lb = input.lb(whichParams);
    ub = input.ub(whichParams);
    
    % create handle to cost function - needs to be instantiated outside
    % parfor loop
    nll_func = @computeNLL;

    % iterate over parameter start points
    % vector to store NNL for each start point
    NLL_sp = inf(1,numSP);
    % vector to store optimized end points for each iteration
    endPoints = nan(numSP,numParams_opt);
    for spi=1:numSP
%     parfor spi=1:numSP
        
        % set initial param start point
        if spi == 1
            X0 = startPoint(whichParams);
        else
            % draw at random between bounds
            X0 = lb + (rand(1,numParams_opt) .* (ub - lb));
        end
        
        % Minimize negative log likelihood of observing data given the
        % parameters
        try
            [X_,NLL] = fmincon(nll_func, X0, [], [], [], [], lb, ub, [], options);
        catch
            % Reselect a random start point if the cost function is
            % undefined at the initial start point
            %originally was && - caused issues when searchSP was false
            %RVT reverting back from || to && to avoid infinite loop
            NLL = nll_func(X0);
            while spi > 1 && ~isfinite(NLL)
                X0 = lb + (rand(1,numParams_opt) .* (ub - lb));
                NLL = nll_func(X0);
            end
            % Rerun minimization
            [X_,NLL] = fmincon(nll_func, X0, [], [], [], [], lb, ub, [], options);
        end
        % Set argmin across initial conditions
        NLL_sp(spi) = NLL;
        endPoints(spi,:) = X_;
    end

    % select best fitting parameters across start point iterations
    [NLL_min,min_idx] = min(NLL_sp);
    endPoint = endPoints(min_idx,:);

    % Store fitted parameters and log likelihood
    output.params(whichParams) = endPoint;
    output.params(~whichParams) = input.startPoint(~whichParams);
    output.logLike = -NLL_min;
       
elseif ~input.doSim
    % Store initial parameters and compute log likelihood
    output.params = startPoint;
    [output.logLike,output.logLikeVec] = computeNLL(startPoint);
    output.logLike = -output.logLike;
    output.logLikeVec = -output.logLikeVec;
else
    % Store initial parameters but don't compute log likelihood
    output.params = startPoint;
end


% Simulate model belief if applicable
if input.doSim

    % Loop through blocks
    output.simBelief = nan(size(input.blocks));
    for b_i=1:numel(endBlockInd)
        
        % Set noise for the given block
        blockType_ind = find(blockTypes==input.blocks(begBlockInd(b_i)));
        if uniqueNoiseBlocks
            noiseStd = output.params(length(output.params)-(numel(blockTypes)-blockType_ind));
        else
            noiseStd = output.params(5);
        end
        blSel=input.blocks==blockTypes(blockType_ind);
        blockSetSize=unique(input.setSizes(blSel)); % Divide buffer size by number of colors!!!
        blockColors=unique(input.colors(blSel)); % Divide buffer size by number of colors!!!

        for c_i=1:length(blockColors)
            sel = blSel & input.colors==blockColors(c_i);
            % Divide buffer size by the task's set size
            [~,selBelief] = bayesWMmodel_forFit(input.outcomes(sel),[],output.params(1)./blockSetSize,output.params(2),output.params(3),noiseStd,output.params(4),input.coordSys);
            output.simBelief(sel) = selBelief;
        end
    end
end


% Function for computing the negative log likelihood of observing all
% blocks of data in the session given the parameters, X
    function [negLogLike,negLogLikeVec] = computeNLL(X)
        uParam(whichParams)=X;
        uParam(~whichParams)=input.startPoint(~whichParams);

        % Loop through blocks: choose noise condition appropriate for block
        negLogLikeVec = nan(size(input.outcomes));
        for bi=1:numel(endBlockInd)
            
            % Set noise for the given block
            blockType_ind = find(blockTypes==input.blocks(begBlockInd(bi)));
            if uniqueNoiseBlocks
                noiseStd = uParam(length(uParam)-(numel(blockTypes)-blockType_ind));
            else
                noiseStd = uParam(5);
            end
            
            blSel=input.blocks==blockTypes(blockType_ind);
            blockSetSize=unique(input.setSizes(blSel)); % Divide buffer size by number of colors!!!
            blockColors=unique(input.colors(blSel)); % Divide buffer size by number of colors!!!
            
            for ci=1:length(blockColors)
                sel = blSel & input.colors==blockColors(ci);
                % Divide buffer size by the task's set size
                [allLikes,~] = bayesWMmodel_forFit(input.outcomes(sel),input.beliefs(sel),uParam(1)./blockSetSize,uParam(2),uParam(3),noiseStd,uParam(4),input.coordSys);
                negLogLikeVec(sel) = -log(allLikes); % Store trial-by-trial -log(likelihood)
            end
        end
        
        % Summate over blocks and convert to negative log likelihood
        negLogLike = sum(negLogLikeVec);
        
        % If applicable, compute and add prior to the negative log
        % likelihood
        if input.usePrior==true && ~isempty(input.paramMean) && ~isempty(input.paramCov)
            %prior = normpdf(uParam(whichParams),input.paramMean(whichParams),input.paramStd(whichParams));
            prior = mvnpdf(uParam(whichParams),input.paramMean(whichParams),input.paramCov(whichParams,whichParams));
            negLogLike = negLogLike - sum(log(prior));
        end
        
        % Notify user if current param settings are unstable 
        if ~isfinite(negLogLike)
            disp('fit_AL_WM_allSetSize warning: reached an infinite negative log likelihood!!')
        end
    end
end

