function plot_paramFitRecovery(predictor,response,paramNames)
% PLOT_PARAMFITRECOVERY(predictor,response,paramNames) Plot m-by-m grid, of
% subplots, each showing a fitted line and R^2 value of the response to 
% predictor variable for the mth dimension.
% 
% predictor: n-by-m matrix (type double) of with n samples (i.e., subjects) 
% and m regressors (i.e., parameters)
%
% response: n-by-m matrix (type double) of with n samples (i.e., subjects) 
% and m response variables (i.e., parameters)


numParam = length(paramNames);

figure
for yi=1:numParam
    for xi=1:numParam
        
        [R,pval] = corr(response(:,yi),predictor(:,xi));
        [B,~,~,~,stats] = regress(response(:,yi),[ones(size(predictor,1),1),predictor(:,xi)]);
        x = predictor(:,xi);
        y = B(2).*x + B(1);
        
        subplot(numParam,numParam,(yi-1)*numParam+xi)
        plot(predictor(:,xi),response(:,yi),'b.')
        hold on
        plot(x,y,'k-')
        hold off
        if yi == numParam
            xlabel(paramNames{xi})
        end
        if xi == 1
            ylabel(paramNames{yi})
        end
        %title({['{\beta}_0=',num2str(B(1)),', {\beta}_1=',num2str(B(2))],['R^2=',num2str(stats(1))]})
        title_str = ['r=',num2str(round(R,2))];
        if pval < 0.05
            title_str = append(title_str,'*');
        end
        title(title_str)
    end
end

set(findall(gcf,'-property','FontSize'),'FontSize',12)

end