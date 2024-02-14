function [output]= fitRedBayesModel_zombies_allSetSize(input)


% Goal -- fit data from all set sizes at once!!!!

% Steps to fit K:
% 1) change  input.whichParams, input.lb, input.ub,input.startPoint
%            to include K -- true, 0, 3, 1.5
% 2) this line "Pm=input.K/input.setSize" should use tParams(indexOfKParam)
%            instead of input.K. 
% 3) make this line 7 elements instead of 6: expTrans=false(1, 6); 
% 4) uncomment these lines: 
 %          like = like * max([1, input.K/input.setSize]) + ...
  %          (1/(2*pi))*(1-max([1, input.K/input.setSize]));
  %         ???? IS IT
  %          SUPPOSED TO BE MIN OR MAX???
% Then run for a single subject, check how it looks to see if things make
% sense... then run for all subjects and create TAC plots.







% input      = structure with following fields:

% PARTS OF DATA:
% PE         = subject prediction errors
% UP         = subject updates
% nn         = noise array (standard deviation of outcomes
% color      = what is the color of zombie type.
% newBlock   = logical marking the first trial of a new block
% setSize   = how many zombie colors were in this block? 


% startPoint  = starting point for fit terms, fixed value for other ones
% whichParams = is a logical array, length=6, specifying which parameters
%               you want to fit.
% usePrior    = do you want to use a set of priors here? they are currently
% hand tuned... so you better make sure that they make sense before setting
% this flag to true.
% input.lb    = lower bound
% input.ub    = upper bound
% input.simFromScratch = if true, model will simulate behavior from
%                         scratch, without conditioning on prev subject prediction history
% input.outcomes = only required if simulating data from scratch:
% input.setSize = set size of this trail


% limits:
%    haz,        lw  ,    stdInt  stdSlope     binUp,  STDofBinDist   % K????
% lb=[ 10e-10 ,  10e-10,      0,      0,            0,            0,  ];
% ub=[ 1-10e-10 ,  1-10e-10,  4,      1,            1,            2,  ];



% when doFit is false, function will be called directly
% without using fmincon to minimize squared error. This should only be used
% only for debugging...


% order for startPoint and whichParams:
%1 = hazard,
%2 = likelihood weight, (1 = optimal, 0 = stupid)
%3 = update standard deviation intercept;                            (LOG*)
%4 = update standard deviation slope (over absolute error magnitude);
%5 = Max prob of binary updates
%6 = Width of binary update distribution
%7 = capacity?

if ~isfield(input, 'skipLike')
    input.skipLike=false(size(input.PE));
end

if ~isfield(input, 'doFit')
    input.doFit=True;
end

% just in case, i guess
wp=logical(input.whichParams);
expTrans=false(1, 7); 


%% Call fmincon to evaluate model likelihood:
%  with a parameters initially set to starting point.

start_point = input.startPoint(wp);
model = @frugFun_zombies_extend;                % this is the function handle to the function that takes the parameters and outputs the thing we want to minimize

output=struct;
if input.doFit
    oldOpts = optimset('fmincon');
    options=optimset(oldOpts, 'maxFunEvals', 100000, 'MaxIter', 100000);
    [estimates, ~, ~, ~, ~, ~, ~] = fmincon(model, start_point, [], [], [], [], input.lb(wp), input.ub(wp), [], options);
    output.params=estimates;
else
    output.params=start_point;
end

[output.negLogLike, output.simUP, output.allLogLike, output.binProb, output.simErrBased_UP] = frugFun_zombies_extend(output.params);


% frugFun accepts some inputs (params) and an array that tells it which
% inputs they are (wp)
    function [negLogLike, simUP, allLogLike, binProb, simErrBased_UP] = frugFun_zombies_extend(params) %input comes from fmincon (initially start_point)
        %Create a set of parameters that are either start point (if we are
        %not fitting them) or current fit parameters if we are:
        tParams=input.startPoint;
        tParams(wp)=params;
        tParams=real(tParams);
        tParams=max([tParams; input.lb]);
        tParams=min([tParams; input.ub]);
        
        % some of the parameters are estimated in log space -- so take the
        % exponent of the parameter in order to do anything:
        tParams(expTrans)=exp(tParams(expTrans));
        
       
        %% LOOP THROUGH SET SIZES and compute LL FOR EACH SET SIZE
        
        setSizes=unique(input.setSize);
        
        totLogLike=0;
        
        for j = 1:length(setSizes) % ONLY ONE SET SIZE? - 3 , should only run once
            sel=input.setSize==setSizes(j); % select data from just one set size
            
            % get all data for just this block:
            justDat=rmfield(input, {'usePrior', 'doFit', 'simFromScratch', 'lb', 'ub', 'startPoint', 'whichParams'});
            justDat=straightStruct(justDat);
            blockDat=selBehav(justDat, sel); % select just data for this block.
            
            blockColors=unique(blockDat.color);
            %keyboard
            
            
            for i = 1:length(blockColors)
                sel2=blockDat.color==blockColors(i);
                colorDat=selBehav(blockDat, sel2); % select just data for this color.
                
                % keyboard
                
                % order for startPoint and whichParams:
                %1 =hazard,
                %2 =likelihood weight,
                %3 =heliVisVar,
                %4 =uncertainty depletion
                %keyboard
                
                NN=colorDat.nn;
                
                
                %   RIGHT NOW THIS WILL ONLY WORK FOR SET SIZE = 1
                [~, ~, ~, errBased_UP]=getTrialVarsFromPEs_cannon( ...
                   NN, colorDat.PE, tParams(1),  colorDat.newBlock, false(size(colorDat.PE)), ...
                    .99, 0, tParams(2), 1, zeros(length(colorDat.newBlock),1), false(length(colorDat.newBlock),1), 2.*pi);
                
                
                % If you make a big error, model will be less certain about exact update.
                
                % Now make variable learning rate in proportion to update:
                % width of update distribution linearly dependent on update size:
                wt=tParams(3)+abs(errBased_UP).*tParams(4); % wt is the standard deviation of a gaussian response distribution
                % ok, we wont let weights go below zero. this is fudgy.
                wt=max([zeros(size(wt))+.0001, wt],[], 2);
                % Likelihood of update is normally distributed around "model"
                % update
                like=normpdf(colorDat.UP-errBased_UP, ...
                    zeros(size(colorDat.UP)), wt);

              
                if length(tParams) == 7
                    like = like * min([1, tParams(7)./unique(colorDat.setSize)]) + ...
                    (1/(2*pi))*(1-min([1, tParams(7)./unique(colorDat.setSize)]));
                end
                
                
                %5 = Max prob of binary updates
                %6 = Width of binary update distribution
                
                
                % probability of no update given NO context error:
                rawP_at_zero=normpdf(0,0,tParams(6));
                pScale = tParams(5)/rawP_at_zero;
                noUpLike=normpdf(errBased_UP, 0, tParams(6)).*pScale; % prob on a normal distribution around 0 update w std tParams(13) scaled by pScale
                

                %  "SOFT" binarization:
                logLike=log( ...
                    (1-noUpLike).*like + ...    % Model predicted update
                    normpdf(colorDat.UP,0,1).*noUpLike);
              
                % THIS HELPS WITH ROBUSTNESS!!! WE are now fitting "average" trials
                % rather than the most weird ones.
                if isfield(input, 'unifMix')
                    % make actual log likelihood a mixture of the one computed
                    % above and a uniform:
                    
                    scaledLL=logLike+log(1-input.unifMix);
                    unifLL=log(input.unifMix)+log(1./300); % random updating.
                    logLike=logsumexp([scaledLL, repmat(unifLL, length(scaledLL), 1)], 2);
                    
                end
                
                % don't let log likelihoods be negative infinity...
                logLike(logLike==-inf)=-10000;
                logLike(logLike==inf)=10000;
                totLogLike=totLogLike+nansum(logLike);
            end
            
        end
        
        % Make total log likelihood negative for minimization function:
        negLogLike=-1.*(totLogLike);

        %ERROR BASED UPDATE ALWAYS WORKS FROM THIS
        if input.simFromScratch
            modPred=nan(size(input.outcomes));
            modPE=nan(size(input.outcomes));
            modUnc=nan(size(input.outcomes));
            modSurprise=nan(size(input.outcomes));

            for ssi = 1:length(setSizes)
                blockSetSize = setSizes(ssi);
                sel_setSize = input.setSize==blockSetSize;
                blockColors = unique(input.color(sel_setSize));
                for ci = 1:length(blockColors)
                    sel = sel_setSize & input.color==blockColors(ci); % select data from just one set size

                    % instantiate block data vectors
                    blockOutcomes = input.outcomes(sel);
                    blockNoise = input.nn(sel);
                    blockPred = nan(sum(sel),1);
                    blockPE = nan(sum(sel),1);
                    blockUnc = nan(sum(sel),1);
                    blockSurprise = nan(sum(sel),1);

                    % initialize dynamic variables
                    blockPred(1) = rand.*2.*pi;
                    blockUnc(1) = .9;

                    % loop through trials and implement learning & choice:
                    for i = 1:length(blockOutcomes)
                        % compute prediction error:
                        blockPE(i)=circ_dist(blockOutcomes(i), blockPred(i)); % subjective prediction error experienced by MODEL
                        
                        % use prediction error to compute model-based variables
                        [modSurp, modUnc, ~, exUP]=getTrialVarsFromPEs_cannon( ...
                            [blockNoise(1),blockNoise(1)], [blockPE(i),blockPE(i)], tParams(1),  [true, false], [false,false], ...
                            blockUnc(i), 0, tParams(2), 1, [0, 0], [false, false], 2.*pi);
                        
                        blockUnc(i+1)=modUnc(2); % choose updated uncertainty, after you see outcome (second index!)
                        blockSurprise(i)=modSurp(1);% choose initial surprise to prediction error, based on previous uncertainty (first index!)
                        modUpdate=exUP(1);
                        
                        % RVT: add variable binary update
                        % probability of no update given NO context error:
                        rawP_at_zero=normpdf(0,0,tParams(6));
                        pScale = tParams(5)/rawP_at_zero;
                        noUpLike=normpdf(modUpdate, 0, tParams(6)).*pScale; % prob on a normal distribution around 0 update w std tParams(13) scaled by pScale
                        
                        Pm=tParams(7)/blockSetSize;
                        if Pm < rand
                            blockPred(i+1)= pi*((rand.*2)-1);
                        elseif rand < noUpLike
                            % don't update belief
                            blockPred(i+1)=blockPred(i);
                        else
                            % update belief
                            % calculate update:
                            wt=tParams(3)+abs(modUpdate).*tParams(4); % wt is the standard deviation of a gaussian response distribution
                            % ok, we wont let weights go below zero. this is fudgy.
                            wt=max([0.0001, wt]);
                            % update:
                            % might have to translate back into range of 0-2 pi
                            blockPred(i+1)=blockPred(i)+normrnd(modUpdate, wt);
                        end
                    end
                    
                    modPred(sel) = blockPred(1:end-1);
                    modPE(sel) = blockPE;
                    modUnc(sel) = blockUnc(1:end-1);
                    modSurprise(sel) = blockSurprise;
                end
            end
            fprintf('modSurprise = %g\n', mean(modUnc))
            output.modPred=circ_dist(modPred, 0); % go from -pi to pi
            output.outcomes=circ_dist(input.outcomes, 0);
            output.modUnc=modUnc;
            output.modSurprise=modSurprise;
            
            simUP=nan(size(input.outcomes));
            binProb=nan(size(input.outcomes));
            % this is for subject errors, not model:
            simErrBased_UP=errBased_UP; %ERROR OCCURING HERE - errBased_UP ALL NANS

            % XXX RVT: not sure what to do with this since ErrBased_UP
            % isn't relevant for normative model belief (right?)
            % ...tried adding something for this above

%             %keyboard
%             % select a random subset of trials for context errors:
%             simErrBased_UP=errBased_UP;
%             
%             % Simulate updates:
%             simUP=simErrBased_UP+normrnd(zeros(size(wt)), wt);
%             
%             % simulate binary updates:
%             
%             % 1) choose trials to binary update
%             % 2) choose whether those trials will be no/total updates
%             
%             % if we did binarize, which way would we go?
%             
%             % create a vector of random numbers:
%             coinFlip=rand(size(noUpLike));
%             isNo=coinFlip<noUpLike;
%             isAll=coinFlip>=noUpLike & coinFlip< noUpLike;
%             fprintf('Mean no update probability = %g, mean tot update prob = %g\n', nanmean(noUpLike), 0)
%             
%             
%             % if we did binarize and we would have gone total, mark it as
%             % such:
%             simUP(isAll)=input.PE(isAll);
%             % if we did binarize and we would have done no update, mark it
%             % as such:
%             simUP(isNo)=0;
%             
%             binProb=[noUpLike, noUpLike];
        end

        logLike=logLike(isfinite(input.UP)&~input.skipLike);
        
        allLogLike=logLike;
%         end
    %    keyboard;
        
        
%         logLike=logLike(isfinite(input.UP)&~input.skipLike);
%         
%         allLogLike=logLike;
        % keyboard
        if ~isfinite(negLogLike)||negLogLike~=real(negLogLike)
            
            keyboard
        end
        
     %   fprintf('negLogLike = %g\n', mean(negLogLike))
    end

end






                  