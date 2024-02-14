function plot_modelCoeff(subjCoeff,modCoeff,modCoeffErr,modLegend)
% Plot subject-by-subject regression coefficients from linear update vs.
% prediction error model (i.e. plots prediction error interactions with CPP
% and RU estimates).
%
% subjCoeff     -- m-by-n matrix (double) of m subjects and n coefficients
% modCoeff      -- k-length cell array of k normative models, each with n
%                  coefficients
% modCoeffErr   -- k-length cell array of n-by-2 (double) of upper and
%                  lower bounds for modCoeff error bars
% modLegend     -- vector of k cells containing the name (string) for each
%                  normative model

getCbColors

numCoeff = size(subjCoeff,2);
numModels = length(modCoeff);
x_subj = repmat(1:numCoeff,size(subjCoeff,1),1) + smartJitter(subjCoeff,0.05,0.06);
% c = [[0.4660 0.6740 0.1880];[0.8500 0.3250 0.0980];[0 0.4470 0.7410];[0.9290 0.6940 0.1250]];

% Plot lines between subject-specific coefficients
figure
plot(x_subj',subjCoeff','color',[0.8 0.8 0.8],'lineWidth',2,'HandleVisibility','off')
hold on

% Plot subject coefficients
for ci=1:numCoeff
    scatter(x_subj(:,ci),subjCoeff(:,ci),[],cbColors(ci,:),'lineWidth',2,'HandleVisibility','off')
end

% Convert to matrix and jitter model coefficients
modCoeff = cell2mat(modCoeff);
x_mod = repmat(1:numCoeff,numModels,1) + smartJitter(modCoeff,0.05,0.06);

% Plot coefficients for each model
for mi=1:numModels
    plot(x_mod(mi,:),modCoeff(mi,:),'-d','markerSize', 9,'markerFaceColor',...
         cbColors(mi,:),'color',cbColors(mi,:),'markerEdgeColor','none')
    
    % Error bars
    if ~isempty(modCoeffErr)
        errorbar(x_mod(mi,:),modCoeff(mi,:),modCoeffErr{mi}(1,:), ...
                 modCoeffErr{mi}(2,:),'lineStyle','none','color', ...
                 cbColors(mi,:),'HandleVisibility','off')
    end
end

line(xlim,[0 0],'color','k','lineStyle',':','HandleVisibility','off')
hold off


ylabel('coefficient')
xlabel('term')
ylim([-2 2])
if ~isempty(modLegend)
    legend(modLegend)
end
xticks(1:numCoeff)
for labi=1:numCoeff
    labels{labi} = ['{\beta}_',num2str(labi)];
end
xticklabels(labels)
set(findall(gcf,'-property','FontSize'),'FontSize',14)

end