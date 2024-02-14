function plot_task(trueMean,outcomes,obsB,outliers,simB,cp,selTrials,coordSys)
% trueMean: vector of m instances of the true expected value, where m is
% the number of trials
% outcomes: vector of m instances of the true outcome
% obsB: vector of m instances of the experimentally observed belief
% outliers: vector of m binary values where 1's denote trials where the agent behaved abnormally
% simB: vector of m instances of the simulated (theoretical) belief
% cp: vector of m instances of binary change point classification
% selTrials: 1x2 vector of integers representing the selected (subject-specific) trials to visualize

% Set default coordinate system if necessary
if nargin<8
    coordSys = 'cartesian';
end

trialNum = 1:length(trueMean);
trialInds = selTrials(1):selTrials(2);
cpInd = find(cp(trialInds)==1); %trial index of change point

% Outlier-corrected observed belief
obsB_corrected = obsB;
outlierInds = find(outliers(trialInds)==1);
for oi=1:numel(outlierInds)
    obsB_corrected(trialInds(outlierInds(oi))) = outcomes(trialInds(outlierInds(oi))-1); % Worst-case sceniorio: subject adopts learning rate of 1
end

% Figure settings
colors = [[0.4940 0.1840 0.5560]; [0.9290 0.6940 0.1250]; [0 0.4 0]];

figure
if isequal(coordSys,'cartesian')
    % Belief trajectories
    plot_names = {'true expected value','trial outcome',"subject's belief","model's belief"};
    subplot(2,1,1)
    plot(trialNum(trialInds),trueMean(trialInds),'lineWidth',2,'color',[0.8 0.8 0.8])
    hold on
    scatter(trialNum(trialInds),outcomes(trialInds),'lineWidth',2,'markerEdgeColor',colors(1,:))
    scatter(trialNum(trialInds),obsB(trialInds),'lineWidth',2,'markerEdgeColor',colors(2,:))
    if ~isempty(outlierInds)
        % Specify outlier trials if they exist
        scatter(trialInds(outlierInds),obsB(trialInds(outliers(trialInds)==1)),'lineWidth',2,...
            'markerEdgeColor',colors(2,:),'markerFaceColor','k')
        plot_names = {'true expected value','trial outcome',"subject's belief",'outlier trial',"model's belief"};
    end
    scatter(trialNum(trialInds),simB(trialInds),'lineWidth',2,'markerEdgeColor',colors(3,:))


    yl = ylim;
    for cp_i=1:numel(cpInd)
        line([trialInds(cpInd(cp_i)) trialInds(cpInd(cp_i))],yl,'lineStyle','--','Color','k')
    end
    hold off
    xlabel('trial number')
    ylabel('location')
    legend(plot_names)
    xlim(selTrials)

    % Distance from true expected value
    subplot(2,1,2)
    plot(trialNum(trialInds),abs(outcomes(trialInds)-trueMean(trialInds)),'color',colors(1,:),'lineWidth',2)
    hold on
    plot(trialNum(trialInds),abs(obsB_corrected(trialInds)-trueMean(trialInds)),'color',colors(2,:),'lineWidth',2)
    plot(trialNum(trialInds),abs(simB(trialInds)-trueMean(trialInds)),'color',colors(3,:),'lineWidth',2)
    yl = ylim;
    for cp_i=1:numel(cpInd)
        line([trialInds(cpInd(cp_i)) trialInds(cpInd(cp_i))],yl,'lineStyle','--','Color','k')
    end
    hold off
    xlim(trialNum(selTrials))
    xlabel('trial number')
    ylabel('dist. from true EV')
    
elseif isequal(coordSys,'polar2D')
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
    trueMean_rad = circ_dist(trueMean(trialInds),0);
    for trial_idx=2:numel(trueMean_rad)
        current_mean_rad = trueMean_rad(trial_idx);
        previous_mean_rad = trueMean_rad(trial_idx-1);
        dist_from_previous = circ_dist(current_mean_rad, previous_mean_rad);
        trueMean_rad(trial_idx) = previous_mean_rad + dist_from_previous;
    end
    interp_locs = trialInds(1):0.001:trialInds(end);
    trueMean_interp = interp1(trialInds,trueMean_rad,interp_locs);
    trueMean_x = cos(trueMean_interp);
    trueMean_y = sin(trueMean_interp);
    trueMean_z = interp_locs;
    
    trialNum_z = trialInds;

    outcomes_dist = circ_dist(trueMean(trialInds),0) - circ_dist(trueMean(trialInds),outcomes(trialInds));
    outcomes_x = cos(outcomes_dist);
    outcomes_y = sin(outcomes_dist);

    obsB_dist = circ_dist(trueMean(trialInds),0) - circ_dist(trueMean(trialInds),obsB(trialInds));
    obsB_x = cos(obsB_dist);
    obsB_y = sin(obsB_dist);

    obsB_outlier_dist = circ_dist(trueMean(trialInds(outlierInds)),0) - circ_dist(trueMean(trialInds(outlierInds)),obsB(trialInds(outlierInds)));
    [obsB_outlier_x, obsB_outlier_y] = pol2cart(obsB_outlier_dist,1);

    simB_dist = circ_dist(trueMean(trialInds),0) - circ_dist(trueMean(trialInds),simB(trialInds));
    [simB_x, simB_y] = pol2cart(simB_dist,1);

    % create cylindrical mesh with radius of 1
    theta = 0:pi/8:2*pi;
    height = trialNum_z(1):round(length(trialNum_z)/5):trialNum_z(end)+1;
    [mesh_theta, mesh_z] = meshgrid(theta, height);

    % Belief trajectories
    plot_names = {'true expected value','trial outcome',"subject's belief","model's belief"};
    subplot(2,1,1);

    % create cylindrical surface
    surf(cos(mesh_theta),sin(mesh_theta),mesh_z,'FaceColor',[0.8 0.8 0.8], ...
         'FaceAlpha',0.1,'EdgeColor',[0.8 0.8 0.8],'HandleVisibility', 'off');
    hold on

    % true mean trajectory
    plot3(trueMean_x,trueMean_y,trueMean_z,'lineWidth',2,'color','k')
    scatter3(outcomes_x,outcomes_y,trialNum_z,'lineWidth',2,'markerEdgeColor',colors(1,:))
    scatter3(obsB_x,obsB_y,trialNum_z,'lineWidth',2,'markerEdgeColor',colors(2,:))
    if ~isempty(outlierInds)
        % Specify outlier trials if they exist
        scatter3(obsB_outlier_x,obsB_outlier_y,trialNum_z,'lineWidth',2,...
            'markerEdgeColor',colors(2,:),'markerFaceColor','k')
        plot_names = {'true expected value','trial outcome',"subject's belief",'outlier trial',"model's belief"};
    end

    % model's belief
    scatter3(simB_x,simB_y,trialNum_z,'lineWidth',2,'markerEdgeColor',colors(3,:))
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
    subplot(2,1,2)
    plot(trialInds,abs(circ_dist(outcomes(trialInds),trueMean(trialInds))),'color',colors(1,:),'lineWidth',2)
    hold on
    plot(trialInds,abs(circ_dist(obsB_corrected(trialInds),trueMean(trialInds))),'color',colors(2,:),'lineWidth',2)
    plot(trialInds,abs(circ_dist(simB(trialInds),trueMean(trialInds))),'color',colors(3,:),'lineWidth',2)
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
end
set(findall(gcf,'-property','FontSize'),'FontSize',14)

end