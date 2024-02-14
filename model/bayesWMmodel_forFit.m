function [B_Like, simulatedData] = bayesWMmodel_forFit(outcomes,B,buffSize,haz,motorStd,noiseStd,cpLikeBias,coordSys)
% Working memory model using change-point discrimination by Bayesian
% estimation over samples in buffer

% Set default coordinate system if necessary
if nargin<7
    coordSys = 'cartesian';
end

% Set number of possible outcomes (specific to experiment)
switch coordSys
    case 'cartesian'
        maxPossOutcome = 300;
    case 'polar'
        maxPossOutcome = 2*pi;
    otherwise
        error("bayesWMmodel_forFit error: didn't recognize mode!! Use mode='cartesian' or 'polar'.")
end
randChoiceProb = 1/(maxPossOutcome); % Uniform pdf over [0,maxPossOutcome]

% Initialize output vectors
B_Like = nan(length(outcomes),1);
simulatedData = nan(length(outcomes),1);

% Simulate data
if isempty(B)
    
    lowBuffChoice = [rand*maxPossOutcome; nan(length(outcomes)-1,1)]; % initial random guess
    highBuffChoice = lowBuffChoice; % the same initial guess 
    
    % Sampling with small buffer size
    minBuff = floor(buffSize);
    for ti=2:length(outcomes)
        if ti<=minBuff
            % Before there are enough samples to fill buffer
            buffer = outcomes(1:ti-1);
        else
            % Randomly assign an extra element to buffer a percentage of
            % the time
            buffer = outcomes(ti-minBuff:ti-1);
        end
        
        % Make a prediction of the current belief given the current buffer
        lowBuffChoice(ti) = makePrediction(buffer,[],maxPossOutcome,randChoiceProb,haz,motorStd,noiseStd,cpLikeBias,coordSys);
    end
    
    % Sampling with large buffer size
    maxBuff = ceil(buffSize);
    for ti=2:length(outcomes)
        if ti<=maxBuff
            % Before there are enough samples to fill buffer
            buffer = outcomes(1:ti-1);
        else
            % Randomly assign an extra element to buffer a percentage of
            % the time
            buffer = outcomes(ti-maxBuff:ti-1);
        end
        
        % Make a prediction of the current belief given the current buffer
        highBuffChoice(ti) = makePrediction(buffer,[],maxPossOutcome,randChoiceProb,haz,motorStd,noiseStd,cpLikeBias,coordSys);
    end
    
    % The final simulated belief is assigned according the rate of high vs
    % low buffer size (as specified by buffSize)
    for ti=1:length(outcomes)
        if rand<ceil(buffSize)-buffSize
            simulatedData(ti) = lowBuffChoice(ti);
        else
            simulatedData(ti) = highBuffChoice(ti);
        end
    end

% Calculate likelihood
else
    
    likeB_min = [randChoiceProb; nan(length(outcomes)-1,1)];
    likeB_max = likeB_min;
    
    % Sampling with small buffer size
    minBuff = floor(buffSize);
    for ti=2:length(outcomes)
        if ti<=minBuff
            % Before there are enough samples to fill buffer
            buffer = outcomes(1:ti-1);
        else
            % Randomly assign an extra element to buffer a percentage of
            % the time
            buffer = outcomes(ti-minBuff:ti-1);
        end
        
        % Make an estimate of the likelihood given the known belief (i.e.,
        % prediction) and current buffer
        [~,likeB_min(ti)] = makePrediction(buffer,B(ti),maxPossOutcome,randChoiceProb,haz,motorStd,noiseStd,cpLikeBias,coordSys);
    end
    
    % Sampling with large buffer size
    maxBuff = ceil(buffSize);
    for ti=2:length(outcomes)
        if ti<=maxBuff
            % Before there are enough samples to fill buffer
            buffer = outcomes(1:ti-1);
        else
            % Randomly assign an extra element to buffer a percentage of
            % the time
            buffer = outcomes(ti-maxBuff:ti-1);
        end
        
        % Make an estimate of the likelihood given the known belief (i.e.,
        % prediction) and current buffer
        [~,likeB_max(ti)] = makePrediction(buffer,B(ti),maxPossOutcome,randChoiceProb,haz,motorStd,noiseStd,cpLikeBias,coordSys);
    end

    % Absolute likelihood is computed as a mixture of likeB_min and 
    % likeB_max max
    maxWt = buffSize-floor(buffSize);
    minWt = 1-maxWt;
    B_Like = likeB_max.*maxWt + likeB_min.*minWt;
    
end
end

    




function [choice,likelihood] = makePrediction(buffer,knownChoice,maxPossOutcome,randChoiceProb,haz,motorStd,noiseStd,cpLikeBias,coordSys)
% Make a new prediction (i.e., choice) or estimate the likelihood of the
% current buffer. If the true choice is known (knownChoice), a likelihood 
% estimate rather than a new choice is made.

    % Set default coordinate system if necessary
    if nargin<8
        coordSys = 'cartesian';
    end

    if isempty(knownChoice)
        likelihood = nan;
        switch length(buffer)
            case 0
                % Uniform random guess
                choice = rand * maxPossOutcome;
            case 1
                % Gaussian random guess centered at previous outcome
                choice = randCartPol(buffer,motorStd,coordSys);
            otherwise
                % Initialize first node
                initMean = buffer(1);
                initProb = randChoiceProb;
                initRunLength = 1;

                % Compute node properties and iterate to next
                for bi=2:length(buffer)
                    [initMean,initProb,initRunLength] = cpPosterior(buffer(bi),randChoiceProb,initMean,initProb,initRunLength,haz,noiseStd,cpLikeBias,coordSys);
                end

                % Gaussian random guess centered at computed expected value
                EV_est = meanCartPol(initMean,1,coordSys,initProb); % Probability-weighted mean of the means
                % Gaussian random guess centered at mode
%                 [~,max_idx] = max(initProb);
%                 EV_est = initMean(max_idx);
                choice = randCartPol(EV_est,motorStd,coordSys); % Random sample
        end
    else
        choice = nan;
        switch length(buffer)
            case 0
                % Prob. given a uniform random guess
                likelihood = randChoiceProb;
            case 1
                % Prob. given a gaussian random guess centered at previous outcome
                likelihood = pdfCartPol(knownChoice,buffer,motorStd,coordSys);
            otherwise
                % Initialize first node
                initMean = buffer(1);
                initProb = randChoiceProb;
                initRunLength = 1;

                % Compute node properties and iterate to next
                for bi=2:length(buffer)
                    [initMean,initProb,initRunLength] = cpPosterior(buffer(bi),randChoiceProb,initMean,initProb,initRunLength,haz,noiseStd,cpLikeBias,coordSys);
                end

                % Prob. given a gaussian random guess centered at computed expected value
                EV_est = meanCartPol(initMean,1,coordSys,initProb); % Probability-weighted mean of the means
                % Prob. given a gaussian random guess centered at mode
%                 [~,max_idx] = max(initProb);
%                 EV_est = initMean(max_idx);
                likelihood = pdfCartPol(knownChoice,EV_est,motorStd,coordSys); 
        end
    end
end






function [newMean,newProb,newRunLength] = cpPosterior(bufferItem,randChoiceProb,initMean,initProb,initRunLength,haz,noiseStd,cpLikeBias,coordSys) 
    % randChoiceProb = the probability of choosing any random value for a
    % trial (from uniform distribution)
    % initMean = should store posterior mean computed on previous node in graph.  
    % initProb = how probable is this actual node to be true, given all
    %               that we know about the data (and other possible nodes)
    % initRunLength  = run length since the most recent changepoint.
    % noiseStd = standard deviation from which outcomes (in the buffer) are
    %               selected
    % haz = hazard rate
    % cpLikeBias = exponent that flattens the probability dist. of data
    %               given the possibility of a CP (or not)
    % coordSys = string identifying whether we are using a cartesian or
    %               polar coordinate system

    
    % Set default coordinate system if necessary
    if nargin<8
        coordSys = 'cartesian';
    end
    
    % Initialize column vectors for new mean, prob, and runLength
    % Add extra element to carry cumulative probablity of a changepoint
    % occuring
    newMean = nan(length(initMean)+1,1);
    newProb = zeros(length(initMean)+1,1);
    newRunLength = nan(length(initMean)+1,1);
    
    for i=1:length(initMean)
        
        % Compute probability given a non-changepoint and set updated mean,
        % probability, and run length for that possibility
        % Sample from normal dist
        nonCPLike = pdfCartPol(bufferItem,initMean(i),sqrt(noiseStd^2/initRunLength(i) + noiseStd^2),coordSys); % Likelihood given a non-changepoint
        nonCPPost = nonCPLike^cpLikeBias * (1-haz); % Posterior probability given a non-changepoint
        
        % Compute joint probability given a changepoint
        CPPost = randChoiceProb^cpLikeBias * haz; % Posterior probability given a changepoint
        
        % Normalize posterior
        nonCPPost = nonCPPost / (nonCPPost+CPPost);
        CPPost = 1-nonCPPost;
        
        %newMean(i) = initMean(i) * initRunLength(i) / (1+initRunLength(i)) + bufferItem * 1/(1+initRunLength(i)); % Running average
        oldMean_newObs = [initMean(i),bufferItem]; % Elements to average
        weight = [initRunLength(i)/(1+initRunLength(i)), 1/(1+initRunLength(i))]; % Weights
        newMean(i) = meanCartPol(oldMean_newObs,2,coordSys,weight); % Running average (computed through a weighted mean update)
        newProb(i) = initProb(i) * nonCPPost; % Total unnormalized posterior prob of a non-changepoint + data
        newRunLength(i) = initRunLength(i) + 1; % Increment run length      

        % Compute joint probability given a changepoint
        % Note that the mean and run length can be updated once since it is
        % the same for all initial expected values
        newProb(end) = newProb(end) + (initProb(i) * CPPost); % Running total unnormalized posterior prob of changepoint + data
    end
    
    % Set updated mean and run length for changepoint condition
    newMean(end) = bufferItem; % Reset expected value
    newRunLength(end) = 1; % Reset run length
    newProb = newProb/sum(newProb); % Normalize
end


function m = meanCartPol(x,dim,mode,w)
% Compute weighted mean of a vector, X, in cartesian or polar coordinates 
% over dimension dim. By default, dim=1 and mode='cartesian'.
    if nargin<4
        if nargin<3
            mode = 'cartesian';
            if nargin<2
                dim = 1;
            end
        end
        w = ones(size(x))/size(x,dim);
    end
    
    switch mode
        case 'cartesian'
            m = sum(w.*x,dim)/sum(w);
        case 'polar'
            m = angle(sum(w.*exp(1i*x),dim)); % get circular weighted mean
        otherwise
            error("bayesWMmodel_forFit::meanCartPol error: didn't recognize mode!! Use mode='cartesian' or 'polar'.")
    end
    
end


function p = pdfCartPol(X,mean,sigma,mode)
% Compute likelihood from the pdf selected by mode. For mode='cartesian', a
% normal distribution is used. For mode='polar', a von Mises distribution 
% is used.

	if nargin<4
        mode = 'cartesian';
    end
    
    switch mode
        case 'cartesian'
            p = normpdf(X,mean,sigma);
        case 'polar'
            p = circ_vmpdf_ov(X,mean,1/sigma^2);
        otherwise
            error("bayesWMmodel_forFit::pdfCartPol error: didn't recognize mode!! Use mode='cartesian' or 'polar'.")
    end
end


function p = randCartPol(mean,sigma,mode)
% Sample randomly from the distribution according to mode. For 
% mode='cartesian', a normal distribution is used. For mode='polar', a von 
% Mises distribution is used.

	if nargin<3
        mode = 'cartesian';
    end
    
    switch mode
        case 'cartesian'
            p = normrnd(mean,sigma);
        case 'polar'
            p = circ_vmrnd(mean,1/sigma^2,1);
        otherwise
            error("bayesWMmodel_forFit::randCartPol error: didn't recognize mode!! Use mode='cartesian' or 'polar'.")
    end
end

