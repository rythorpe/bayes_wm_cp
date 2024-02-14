function plot_task_multicolor(trialColor,trueMean,outcomes,obsB,outliers,simB,cp,selTrials,coordSys)
% similar to plot_task.m, but for many sets (each with its own color) at
% the same time; 'cartesian' coordinate system isn't an option here;
% only works properly with coorSys='polar3D'
%
% trialcolor: vector of m instances of the color group associated with each
%             trial, where m is the number of trials
% trueMean: vector of m instances of the true expected value
% outcomes: vector of m instances of the true outcome
% obsB: vector of m instances of the experimentally observed belief
% outliers: vector of m binary values where 1's denote trials where the agent behaved abnormally
% simB: vector of m instances of the simulated (theoretical) belief
% cp: vector of m instances of binary change point classification
% selTrials: 1x2 vector of integers representing the selected (subject-specific) trials to visualize

% Set default coordinate system if necessary
if nargin<8
    coordSys = 'polar3D';
end

% import color palette
getCbColors;

% get trial metadata
allTrialInds = 1:length(trialColor);
trialMask = zeros([length(trialColor),1]);
trialInds = selTrials(1):selTrials(2);
trialMask(trialInds) = 1;
cpInd = find(cp(trialInds)==1); %trial index of change point
colorLabels = unique(trialColor);
n_colors = numel(colorLabels);

% Outlier-corrected observed belief
obsB_corrected = obsB;
outlierInds = find(outliers(trialInds)==1);
for oi=1:numel(outlierInds)
    obsB_corrected(trialInds(outlierInds(oi))) = outcomes(trialInds(outlierInds(oi))-1); % Worst-case sceniorio: subject adopts learning rate of 1
end

% Figure settings
% colors = [[0.4940 0.1840 0.5560]; [0.9290 0.6940 0.1250]; [0 0.4 0]];
colors = cbColors;

figure
if isequal(coordSys,'polar2D')
    trueMean_rad = circ_dist(trueMean(trialInds),0);
    outcomes_dist = circ_dist(outcomes(trialInds),trueMean_rad);
    obsB_dist = circ_dist(obsB(trialInds),trueMean_rad);
    simB_dist = circ_dist(simB(trialInds),trueMean_rad);

    % trial trajectories
    plot_names = {'true expected value','trial outcome',"subject's belief","model's belief"};
    subplot(2,1,1)
    % true mean
    plot(trialInds,trueMean_rad,'lineWidth',2,'color',[0.8 0.8 0.8])
    hold on
    % outcome
    scatter(trialInds,trueMean_rad+outcomes_dist,'lineWidth',2,'markerEdgeColor',colors(1,:))
    % observer's belief
    scatter(trialInds,trueMean_rad+obsB_dist,'lineWidth',2,'markerEdgeColor',colors(2,:))
    if ~isempty(outlierInds)
        % Specify outlier trials if they exist
        obsB_outlier_dist = circ_dist(obsB(trialInds(outlierInds)),trueMean_rad);
        scatter(trialInds(outlierInds),trueMean_rad+obsB_outlier_dist,'lineWidth',2,...
            'markerEdgeColor',colors(2,:),'markerFaceColor','k')
        plot_names = {'true expected value','trial outcome',"subject's belief",'outlier trial',"model's belief"};
    end
    % model's belief
    scatter(trialInds,trueMean_rad+simB_dist,'lineWidth',2,'markerEdgeColor',colors(3,:))
    
    % change points
    yl = ylim;
    for cp_i=1:numel(cpInd)
        line([trialInds(cpInd(cp_i)) trialInds(cpInd(cp_i))],yl,'lineStyle','--','Color','k')
    end
    hold off
    xlim(selTrials)
    xlabel('trial number')
    ylim([-3*pi/2,3*pi/2])
    yticks([-pi,0,pi])
    yticklabels({'-\pi','0','\pi'})
    ylabel('location')
    legend(plot_names)
    
    % Distance from true expected value
    subplot(2,1,2)
    plot(trialInds,abs(circ_dist(outcomes(trialInds),trueMean_rad)),'color',colors(1,:),'lineWidth',2)
    hold on
    plot(trialInds,abs(circ_dist(obsB_corrected(trialInds),trueMean_rad)),'color',colors(2,:),'lineWidth',2)
    plot(trialInds,abs(circ_dist(simB(trialInds),trueMean_rad)),'color',colors(3,:),'lineWidth',2)
    yl = ylim;
    for cp_i=1:numel(cpInd)
        line([trialInds(cpInd(cp_i)) trialInds(cpInd(cp_i))],yl,'lineStyle','--','Color','k')
    end
    hold off
    xlim(selTrials)
    xlabel('trial number')
    ylim([0,pi])
    yticks([0,pi])
    yticklabels({'0','\pi'})
    ylabel('dist. from true EV')

elseif isequal(coordSys,'polar3D')
    
    % create cylindrical mesh with radius of 1
    theta = 0:pi/8:2*pi;
    height = trialInds(1):round(length(trialInds)/5):trialInds(end)+1;
    [mesh_theta, mesh_z] = meshgrid(theta, height);

    % plot cylindrical surface
    subplot(1,2,1);
    surf(cos(mesh_theta),sin(mesh_theta),mesh_z,'FaceColor',[0.8 0.8 0.8], ...
             'FaceAlpha',0.3,'EdgeColor',[0.8 0.8 0.8],'HandleVisibility', 'off');
    
    % plot stuff on top of surface
    for color_idx=1:n_colors
        colorLabel = colorLabels(color_idx);
        trial_sel = (trialColor == colorLabel) & trialMask;
        if sum(trial_sel) < 1
            continue
        end
        color = colors(color_idx,:);
        trialNum_z = allTrialInds(trial_sel);
        
        trueMean_rad = circ_dist(trueMean(trial_sel),0);
        for trial_idx=2:numel(trueMean_rad)
            current_mean_rad = trueMean_rad(trial_idx);
            previous_mean_rad = trueMean_rad(trial_idx-1);
            dist_from_previous = circ_dist(current_mean_rad, previous_mean_rad);
            trueMean_rad(trial_idx) = previous_mean_rad + dist_from_previous;
        end
        interp_locs = trialInds(1):0.001:trialInds(end);
        trueMean_interp = interp1(trialNum_z,trueMean_rad,interp_locs);
        trueMean_x = cos(trueMean_interp);
        trueMean_y = sin(trueMean_interp);
        trueMean_z = interp_locs;
    
        outcomes_dist = circ_dist(trueMean(trial_sel),0) - circ_dist(trueMean(trial_sel),outcomes(trial_sel));
        outcomes_x = cos(outcomes_dist);
        outcomes_y = sin(outcomes_dist);
    
        obsB_dist = circ_dist(trueMean(trial_sel),0) - circ_dist(trueMean(trial_sel),obsB(trial_sel));
        obsB_x = cos(obsB_dist);
        obsB_y = sin(obsB_dist);
        
        outlier_sel = trial_sel & outliers;
        obsB_outlier_dist = circ_dist(trueMean(outlier_sel),0) - circ_dist(trueMean(outlier_sel),obsB(outlier_sel));
        [obsB_outlier_x, obsB_outlier_y] = pol2cart(obsB_outlier_dist,1);
    
        simB_dist = circ_dist(trueMean(trial_sel),0) - circ_dist(trueMean(trial_sel),simB(trial_sel));
        [simB_x, simB_y] = pol2cart(simB_dist,1);
    
        plot_names = {'true expected value','trial outcome',"subject's belief","model's belief"};
        subplot(1,2,1);
        hold on

        % true mean trajectory
        plot3(trueMean_x,trueMean_y,trueMean_z,'lineWidth',2,'color',color)

        % outcome
        scatter3(outcomes_x,outcomes_y,trialNum_z,'lineWidth',2,'markerEdgeColor',color,'marker','o')

        % observer's belief
        scatter3(obsB_x,obsB_y,trialNum_z,'lineWidth',2,'markerEdgeColor',color,'marker','.')
        if sum(outlier_sel) > 0
            % Specify outlier trials if they exist
            scatter3(obsB_outlier_x,obsB_outlier_y,trialInds(outlier_sel), ...
                     'lineWidth',2,'markerEdgeColor',color,'markerFaceColor','k')
            plot_names = {'true expected value','trial outcome',"subject's belief",'outlier trial',"model's belief"};
        end
    
        % model's belief
        scatter3(simB_x,simB_y,trialNum_z,'lineWidth',2,'markerEdgeColor',color,'marker','x')
        hold off

        zlim(selTrials)
        zlabel('trial number')
        xlim([-1,1])
        xticks(0)
        xticklabels({})
        ylim([-1,1])
        yticks(0)
        yticklabels({})
        legend(plot_names)
        
        % Distance from true expected value
        subplot(1,2,2)
        hold on
        % dist_from_mean_outcome = abs(circ_dist(outcomes(trial_sel),trueMean(trial_sel)));
        % set trials belonging to other colors to zero for width consistency in bar plot
        dist_from_mean_outcome = zeros(length(allTrialInds), 1);
        dist_from_mean_outcome(trial_sel) = abs(circ_dist(outcomes(trial_sel),trueMean(trial_sel)));
%         stairs(trialNum_z,dist_from_mean_outcome,'Color',color,'lineWidth',2,'lineStyle','-')
%         bar(trialNum_z,dist_from_mean_outcome,'barWidth',1/n_colors,'faceColor',color,'faceAlpha',0.5,'edgeColor','none')
        bar(allTrialInds,dist_from_mean_outcome,'barWidth',1,'faceColor',color,'faceAlpha',0.5,'edgeColor','none')
        dist_from_mean_obs = abs(circ_dist(obsB_corrected(trial_sel),trueMean(trial_sel)));
        plot(trialNum_z,dist_from_mean_obs,'lineStyle','none','Color',color,'marker','.')
        dist_from_mean_sim = abs(circ_dist(simB(trial_sel),trueMean(trial_sel)));
        plot(trialNum_z,dist_from_mean_sim,'lineStyle','none','Color',color,'marker','x')
        ylim([0,pi])

        % change points
        color_cp_idxs = find(cp & trial_sel);
        for cp_i=1:numel(color_cp_idxs)
            cp_idx = color_cp_idxs(cp_i);
            line([allTrialInds(cp_idx),allTrialInds(cp_idx)],ylim, ...
                 'Color',color,'lineWidth',1,'lineStyle',':')
        end
        hold off
        xlim(selTrials)
        xlabel('trial number')
        yticks([0,pi])
        yticklabels({'0','\pi'})
        ylabel('dist. from true EV')
    end
    set(gca,'view',[90 -90]) % rotate about y=x
end
set(findall(gcf,'-property','FontSize'),'FontSize',14)

end