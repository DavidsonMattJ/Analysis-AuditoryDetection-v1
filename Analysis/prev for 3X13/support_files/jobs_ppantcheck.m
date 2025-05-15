%% This script is called form JOBS_import_preprocess and plots alpha and beta values of ppants to show any anomolies

set_myOnlineDirectories;

open("allppantsTable.mat")

open("allppantsData");

ppants = allppantsData.ppant;

% Plot Alpha values across all speeds

alphasAll = allppantsData.alphaAll;
alphasSlow = allppantsData.alphaSlow;
alphasNorm = allppantsData.alphaNorm;

figure(1); clf;
set(gcf,'units','normalized','position', [0 0 2000 1500]);

% subplot ('position', [100, 100, 2000, 1000]);

plot(ppants, alphasAll, 'o-k');
hold on
plot(ppants, alphasSlow, 'o-b');
plot(ppants, alphasNorm, 'o-r');

title("Threshold Values, per Participant");
xlabel("Participant")
ylabel("Threshold of Detection (dB relative to noise level)")
%ylim([min(alphasSlow), max(alphasNorm)]);

cd(figdir);
cd('group_fitresults');
print('-dpng', 'perparticipant_alphavalues');
%%
% Plot Beta Values 

betasAll = allppantsData.betaAll;
betasSlow = allppantsData.betaSlow;
betasNorm = allppantsData.betaNorm;

figure(1);
clf
plot(ppants, betasAll, 'o-k');
hold on
plot(ppants, betasSlow, 'o-b');
plot(ppants, betasNorm, 'o-r');

title("Beta Values, per participant");
xlabel("Participant")
ylabel("Beta")


cd(figdir);
cd('group_fitresults');
print('-dpng', 'perparticipant_betavalues');
%ylim([0, max(betasNorm)]);





%% Group Means

% groupMeansTable = NaN(2,3);
% groupMeansTable = array2table(groupMeansTable, 'VariableNames', {'MeanAll', 'MeanSlow', 'MeanNorm'}, 'RowNames', { 'Alpha', 'Beta'});
% 
% 
% grpMeanAlphaAll = mean(allppantsData.alphaAll);
% groupMeansTable.MeanAll(1) = grpMeanAlphaAll;
% 
% grpMeanAlphaSlow = mean(allppantsData.alphaSlow);
% groupMeansTable.MeanSlow(1) = grpMeanAlphaSlow;
% 
% grpMeanAlphaNorm = mean(allppantsData.alphaNorm);
% groupMeansTable.MeanNorm(1) = grpMeanAlphaNorm;
% 
% grpMeanBetaAll = mean(allppantsData.betaAll);
% groupMeansTable.MeanAll(2) = grpMeanBetaAll;
% 
% grpMeanBetaSlow = mean(allppantsData.betaSlow);
% groupMeansTable.MeanSlow(2) = grpMeanBetaSlow;
% 
% grpMeanBetaNorm = mean(allppantsData.betaNorm);
% groupMeansTable.MeanNorm(2) = grpMeanBetaNorm;


%% Calculate Group Means

groupmeanAlphas = [mean(alphasAll), mean(alphasSlow), mean(alphasNorm)];

groupmeanBetas = [mean(betasAll), mean(betasSlow), mean(betasNorm)];

SEMs = [std(alphasAll)/sqrt(height(alphasAll)) , std(alphasSlow)/sqrt(height(alphasSlow)) , std(alphasNorm)/sqrt(height(alphasNorm))];

GroupStats = [groupmeanAlphas; groupmeanBetas; SEMs]

%% Plot

conditions = [ 0 , 1 , 2];
labels = {'All', 'Slow', 'Norm'};

figure(1); clf;
bar(groupmeanAlphas);
hold on
errorbar(groupmeanAlphas,SEMs);

ylabel("Threshold (dB)")

title("Group Mean Detection Threshold, by Walk Speed")

hold off

cd(figdir);
cd('group_fitresults');
print('-dpng', 'Group_alphavalues');



%% Threshold change

allppantsData.deltaThresh = allppantsData.alphaNorm - allppantsData.alphaSlow;

% Plot Threshold Change

figure(1);
clf
plot(ppants, allppantsData.deltaThresh);

title("Detection Threshold Change, Per participant");
xlabel("Participant");
ylabel("Change in Threshold (dB) when forced to walk slowly");

cd(figdir);
cd('group_fitresults');
print('-dpng', 'perparticipant_alphaChange');



% %% Mean Threshold Change
% 
% meanDeltaThresh = mean(allppantsData.deltaThresh);
% 
% figure
% plot(ppants, meanDeltaThresh)




