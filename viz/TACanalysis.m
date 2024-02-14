function [output]= TACanalysis(noCorrect,Polar,goodSub)
%CREATES PLOTS ANALYSING TRIAL AFTER CHANGEPOINT INFO
%1- TAC vs Learning Rate for set size 1,2 and 3
%2- TAC vs Guessing for set size 1, 2, and 3
%GOOD SUB determines which row you pull out of params - which subjects has
%good data

Params1=reshape(cell2mat(noCorrect),[],5,3);

Params2=reshape(cell2mat(Polar),[],5,3);

for i=1:size(Params2,1)
    for j=1:size(Params2,2)
        for k=1:size(Params2,3)
            Params(i,j,k)=(Params2(i,j,k)+Params1(i,j,k))./2;
        end
    end
end
%Params=reshape(cell2mat(regress_params),[],5,3);

Params_1=Params(:,:,1);
Params_1=Params_1(goodSub,2:5);
setSizeLabels={'1', '2','3','4'};
xJit=normrnd(0, .1, size(Params_1, 1), 1);
getCbColors

%PLOTTING TAC vs LR for set size 1,2, and 3
figure
subplot(2,3,1)
hold on
for i = 1:size(Params_1, 2)
    plot(i+xJit,Params_1(:,i),'o','markerSize', 5, 'markerFaceColor', cbColors(i+1,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
    ste=std(Params_1(:,i))./sqrt(length(Params_1(:,i)));
    plot([i, i], [mean(Params_1(:,i))+ste, mean(Params_1(:,i))-ste], '-k')  
    plot(i, mean(Params_1(:,i)), 'd', 'markerSize', 15, 'markerFaceColor', cbColors(i+1,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
end
set(gca, 'xtick', 1:4, 'xticklabel', setSizeLabels, 'box', 'off')
ylim([0 1.2])
title('Set Size 1')
xlabel('TAC')
ylabel('Learning Rate')

Params_2=Params(:,:,2);
Params_2=Params_2(goodSub,2:5);
xJit=normrnd(0, .1, size(Params_2, 1), 1);
getCbColors

subplot(2,3,2)
hold on
for i = 1:size(Params_2, 2)
    plot(i+xJit,Params_2(:,i),'o','markerSize', 5, 'markerFaceColor', cbColors(i+1,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
    ste=std(Params_2(:,i))./sqrt(length(Params_2(:,i)));
    plot([i, i], [mean(Params_2(:,i))+ste, mean(Params_2(:,i))-ste], '-k')  
    plot(i, mean(Params_2(:,i)), 'd', 'markerSize', 15, 'markerFaceColor', cbColors(i+1,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
end
set(gca, 'xtick', 1:4, 'xticklabel', setSizeLabels, 'box', 'off')
ylim([0 1.2])
title('Set Size 2')
xlabel('TAC')
ylabel('Learning Rate')
Params_3=Params(:,:,3);
Params_3=Params_3(goodSub,2:5);

xJit=normrnd(0, .1, size(Params_3, 1), 1);
getCbColors

subplot(2,3,3)
hold on
for i = 1:size(Params_3, 2)
    plot(i+xJit,Params_3(:,i),'o','markerSize', 5, 'markerFaceColor', cbColors(i+1,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
    ste=std(Params_3(:,i))./sqrt(length(Params_3(:,i)));
    plot([i, i], [mean(Params_3(:,i))+ste, mean(Params_3(:,i))-ste], '-k')  
    plot(i, mean(Params_3(:,i)), 'd', 'markerSize', 15, 'markerFaceColor', cbColors(i+1,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
end
set(gca, 'xtick', 1:4, 'xticklabel', setSizeLabels, 'box', 'off')
title('Set Size 3')
ylabel('Learning Rate')
xlabel('TAC')
ylim([0 1.2])

 for i=1:size(Params2,1)
    for j=1:size(Params2,2)
        for k=1:size(Params2,3)
            Params(i,j,k)=(Params2(i,j,k)-Params1(i,j,k));
        end
    end
 end

Params_1=Params(:,:,1);
Params_1=Params_1(goodSub,2:5);

xJit=normrnd(0, .1, size(Params_1, 1), 1);
getCbColors

%Plotting TAC vs Guessing for set size 1, 2, and 3
%Set size 1
subplot(2,3,4)
hold on
for i = 1:size(Params_1, 2)
    plot(i+xJit,Params_1(:,i),'o','markerSize', 5, 'markerFaceColor', cbColors(i+1,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
    ste=std(Params_1(:,i))./sqrt(length(Params_1(:,i)));
    plot([i, i], [mean(Params_1(:,i))+ste, mean(Params_1(:,i))-ste], '-k')  
    plot(i, mean(Params_1(:,i)), 'd', 'markerSize', 15, 'markerFaceColor', cbColors(i+1,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
end
set(gca, 'xtick', 1:4, 'xticklabel', setSizeLabels, 'box', 'off')
ylim([-0.05 1.2])
xlabel('TAC')
ylabel('Guessing')

%Set size 2
Params_2=Params(:,:,2);
Params_2=Params_2(goodSub,2:5);

xJit=normrnd(0, .1, size(Params_2, 1), 1);
getCbColors

subplot(2,3,5)
hold on
for i = 1:size(Params_2, 2)
    plot(i+xJit,Params_2(:,i),'o','markerSize', 5, 'markerFaceColor', cbColors(i+1,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
    ste=std(Params_2(:,i))./sqrt(length(Params_2(:,i)));
    plot([i, i], [mean(Params_2(:,i))+ste, mean(Params_2(:,i))-ste], '-k')  
    plot(i, mean(Params_2(:,i)), 'd', 'markerSize', 15, 'markerFaceColor', cbColors(i+1,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
end
set(gca, 'xtick', 1:4, 'xticklabel', setSizeLabels, 'box', 'off')
ylim([-0.05 1.2])
xlabel('TAC')
ylabel('Guessing')

%Set size 3
Params_3=Params(:,:,3);
Params_3=Params_3(goodSub,2:5);

xJit=normrnd(0, .1, size(Params_3, 1), 1);
getCbColors

subplot(2,3,6)
hold on
for i = 1:size(Params_3, 2)
    plot(i+xJit,Params_3(:,i),'o','markerSize', 5, 'markerFaceColor', cbColors(i+1,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
    ste=std(Params_3(:,i))./sqrt(length(Params_3(:,i)));
    plot([i, i], [mean(Params_3(:,i))+ste, mean(Params_3(:,i))-ste], '-k')  
    plot(i, mean(Params_3(:,i)), 'd', 'markerSize', 15, 'markerFaceColor', cbColors(i+1,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
end
set(gca, 'xtick', 1:4, 'xticklabel', setSizeLabels, 'box', 'off')
ylim([-0.05 1.2])
xlabel('TAC')
ylabel('Guessing')


end