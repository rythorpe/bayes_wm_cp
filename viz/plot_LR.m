function plot_LR(outcomes, predictions, blkBeg)
% plot learning rate as a function of prediction error in polar coordinates
% outcomes: m-by-n matrix of trial outcomes, with m trials and n different
%           conditions (each of which will be computed separately and
%           plotted as a different trace)
% predictions: m-by-n matrix of trial predictions
% blkBeg: m-by-n matrix of logicals denoting which trials are the first of
%         a given block

n_bins = 15;
bin_width = pi/n_bins;
PE_bins = 0:bin_width:pi;
n_conds = size(outcomes,2);

figure
for ci=1:n_conds
    [~,UP_noCorrect,PE] = computeLR(outcomes(:,ci),predictions(:,ci), ...
                                    blkBeg(:,ci),'polarNoCorrect');
    [~,UP_polar,~] = computeLR(outcomes(:,ci),predictions(:,ci), ...
                               blkBeg(:,ci),'polar');

    [~, sort_idxs] = sort(abs(PE));
    PE = PE(sort_idxs);
    UP_noCorrect = UP_noCorrect(sort_idxs);
    UP_polar = UP_polar(sort_idxs);
    
    LR = nan(1,n_bins);
    for bi=1:n_bins
        sel_bin = abs(PE) >= PE_bins(bi) & abs(PE) < PE_bins(bi+1);

        Y = UP_noCorrect(sel_bin);
        X = [ones(size(Y)),PE(sel_bin)];
        B = regress(Y,X);
        LR_noCorrect = B(2);

        Y = UP_polar(sel_bin);
        B = regress(Y,X);
        LR_polar = B(2);

        % note: using the average of the polar and non-corrected polar LR
        % has an interesting effect that helps us visualize what's going on
        LR(bi) = (LR_noCorrect + LR_polar) ./ 2;
%         LR(bi) = LR_polar;
%         LR(bi) = LR_noCorrect;
    end

    color = [1, 1, 1] .* 0.9 .* (ci/n_conds);
    hold on
    PE_binCenters = PE_bins(1:end-1) + bin_width ./ 2;
    plot(PE_binCenters,LR,lineWidth=2,color=color)
    hold off
end

xlim([0,pi])
xticks([0,pi])
xticklabels({'0','\pi'})
ylim([0,1.02])
yticks([0,0.5,1])
xlabel('absolute error magnitude')
ylabel('learning rate')
set(findall(gcf,'-property','FontSize'),'FontSize',14)

end