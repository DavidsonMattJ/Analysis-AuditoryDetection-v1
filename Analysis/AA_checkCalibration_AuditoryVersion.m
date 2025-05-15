%quick calibration check:

%quickly load raw csv (summary), and plot the results of calibration.

%%
%%%%%% Auditory& EEG version
set_myOnlineDirectories_AUD

% note now updated for the split dataset (with and without eye).
sancheckplots=0; % toggle to show FA per trial etc.
cd(datadir)

%%
GFX=[]; % stash the DVs for a quick super plot at the end.

%%  Import from csv. summary data.

%frame by frame first:
%PC
pfols = dir([pwd filesep '*trialsummary.csv']);
nsubs= length(pfols);
tr= table([1:length(pfols)]',{pfols(:).name}' );
disp(tr)

%
%% Per csv file, import and wrangle into Matlab Structures, and data matrices:
for ippant =1:length(pfols)
  
cd(datadir)
% cd ../
pfols = dir([pwd filesep '*trialsummary.csv']);
nsubs= length(pfols);
%
filename = pfols(ippant).name;
  
    %extract name&date from filename:
    ftmp = find(filename =='_');
    subjID = filename(1:ftmp(end)-1);
    %read table
    opts = detectImportOptions(filename,'NumHeaderLines',0);
    T = readtable(filename,opts);
    rawSummary_table = T;
    disp(['Preparing participant ' T.participant{1} ]);
    %%
    
    savename = [subjID '_summary_data'];
    % summarise relevant data:
   %  % calibration is performed after every dual target presented
   % practIndex = find(T.isPrac ==1);
   exprows = find(T.isPrac==0);
   exprows= find(T.block>0);
   %  nPracBlocks = unique(T.block(practIndex));
   %  nPrac=length(nPracBlocks);
   %  ss= size(practIndex,1);
 
    %%
    % disp([subjID ' has ' num2str((T.trial(practIndex(end)) +1)) ' practice trials']);
    %extract the rows in our table, with relevant data for assessing
    %calibration:
   
     %% figure 1
    figure(1);  clf; 
    set(gcf, 'color', 'w', 'units', 'normalized', 'position', [0.1 0.1 .8 .8]);

    speedCols={'b', [1, 171/255, 64/255]}; % 
    hold on;
   nBlocks = unique(T.block);
   for iblock = 1:length(nBlocks)
       usetrials = find(T.block==nBlocks(iblock));

       %calculate accuracy:
       calibData = T.correctResponse(usetrials);
       calibAcc = zeros(1, length(calibData));
       for itarg=1:length(calibData)
           tmpD = calibData(1:itarg);
           calibAcc(itarg) = nansum(tmpD)/length(tmpD);
       end

       %retain contrast values:
       calibTone= T.toneAmp(usetrials);

        calibTone(calibTone==-100)=nan;
       %
       subplot(2,10,iblock);
       if iblock==1
        title([subjID ]);
         hold on;
       end
       plot(calibTone, 'o-', 'color', [.4 .4 .4]);
       hold on; ylabel('toneAmp')
       xlabel('target count');
       ylim([-100 0]);
       % add the final contr values:
       %  expContrasts= unique(T.targContrast(exprows));
       % for ic= 1:length(expContrasts)
       % plot(xlim, [expContrasts(ic), expContrasts(ic)], ['r:']);
       % end
       hold on;
       subplot(2,10,iblock+10)
%        yyaxis right
if iblock==1
    col = 'k';
else
    col=speedCols{T.walkSpeed(usetrials(1))};
end
        plot(calibAcc, 'o-', 'color', col);
       ylabel('Accuracy');
       hold on; 
       ylim([0 1])
       xlabel('target count');
       title(['Block ' num2str(iblock)])
shg
   end
   
   
   
   %%
   try cd([figdir filesep 'Calibration']);
   catch mkdir([figdir filesep 'Calibration']);
        cd([figdir filesep 'Calibration']);
   end
   print('-dpng', [subjID ' quick summary '])
    %%
    
    %%
%     legend(lg, 'autoupdate', 'off');
    
  %% now plot accuracy at diff speeds.
    %find standing and walking trials.
    wkSlow= find(T.walkSpeed==1);
    wkNorm= find(T.walkSpeed==2);
   
    %intersect of experiment, and relevant condition:
    walkSrows = intersect(wkSlow, exprows);
    walkNrows = intersect(wkNorm, exprows);
    %% accuracy:
    walkAccSlow = sum(T.correctResponse(walkSrows))/ length(walkSrows);    
    walkAccNorm = sum(T.correctResponse(walkNrows))/ length(walkNrows);
% criterion and sensitivity. 

% we need a row for H, M, FA, CR.
%% add a column for SDT categories
T.SDTcat = nan(size(T,1),1);
sigPresent = find(T.tonePresent);
sigAbsent= find(T.tonePresent==0);
respPres = find(T.targResponse==1);
respAbs = find(T.targResponse==0);
respMissing= find(T.targResponse==-10);

%fill table.
%True positive (Hits)
TP = intersect(sigPresent,respPres);
%True Negative (CR)
TN = intersect(sigAbsent,respAbs);
% Miss
Misses = intersect(sigPresent, respAbs); 
FA = intersect(sigAbsent, respPres);
%1,2,3,4. H,M,FA,CR
T.SDTcat(TP)=1;
T.SDTcat(Misses)=2;
T.SDTcat(FA)=3;
T.SDTcat(TN)=4;

%HR and FAr per type:
SDT_slow= T.SDTcat(walkSrows);
SDT_norm= T.SDTcat(walkNrows);
Pres_slow = nansum(SDT_slow==1) + nansum(SDT_slow==2);
Pres_norm= nansum(SDT_norm==1) + nansum(SDT_norm==2);
Neg_slow = nansum(SDT_slow==3) + nansum(SDT_slow==4);
Neg_norm= nansum(SDT_norm==3) + nansum(SDT_norm==4);

%hit rate
HR_N= nansum(SDT_norm==1)/Pres_norm;
HR_S= nansum(SDT_slow==1)/Pres_slow;

%falsealarm rate
FAR_N= nansum(SDT_norm==3)/Neg_norm;
FAR_S= nansum(SDT_slow==3)/Neg_slow;


vars = [HR_N, HR_S, FAR_N, FAR_S];
vars(vars==1)=.999; % avoid ceiling
vars(vars==0)=.001; % avoid floor
[HR_N, HR_S, FAR_N, FAR_S] = deal(vars(1), vars(2),vars(3),vars(4));

%dprime:
dp_N = norminv(HR_N)- norminv(FAR_N);
dp_S = norminv(HR_S)- norminv(FAR_S);

%criterion
crit_N = -0.5*(norminv(HR_N)+ norminv(FAR_N));
crit_S = -0.5*(norminv(HR_S)+ norminv(FAR_S));

%%
       
    figure(2); clf;
    set(gcf, 'color', 'w', 'units', 'normalized', 'position', [0.1 0.1 .8 .8]);
    
    
   
    subplot(231)
    bar([walkAccSlow, walkAccNorm])
    hold on;  
    title([subjID ]);
    set(gca,'xticklabels', {'slow', 'normal'})
    ylabel('meanAccuracy');
    xlabel('Walk Speed')
    ylim([.2 1]);
    hold on;
        
    text(1-.25, walkAccSlow*1.05, [num2str(round(walkAccSlow,2))]);
    text(2-.25, walkAccNorm*1.05, [num2str(round(walkAccNorm,2))]);
    
    subplot(232);
     bar([HR_S, HR_N; FAR_S, FAR_N]);
    hold on;  
    title(subjID);
    set(gca,'xticklabels', {'Hit Rate', 'False alarm rate'})
    ylabel('Rate');
    xlabel('DV')
    ylim([0 1]);
    legend('Slow', 'Normal');
    hold on;

    %stash for super plot (GFX at end).
    GFX(ippant).slow_acc= walkAccSlow;
    GFX(ippant).norm_acc= walkAccNorm;
    GFX(ippant).slow_HR= HR_S;
    GFX(ippant).norm_HR= HR_N;
    GFX(ippant).slow_FAr= FAR_S;
    GFX(ippant).norm_FAr = FAR_N;
        %%
    %% also RT  
        subplot(234)
  
    allRTs= T.targRT - T.targOnset;
    %remove negatives(these were no resp).
    allRTs(allRTs<0) = nan;
    
     % walkSRT= nansum(allRTs(walkSrows))/ length(walkSrows);   
     % walkRT = nansum(allRTs(walkNrows))/ length(walkNrows);
     walkSRT = nanmean(allRTs(walkSrows));
     walkRT = nanmean(allRTs(walkNrows));
  
    bar([walkSRT, walkRT], 'Facecolor', 'r')
    shg
    title(subjID);
    set(gca,'xticklabels', {'slow', 'normal'})
    ylabel('meanRT');
    
    text(1-.25, walkSRT*1.05, [num2str(round(walkSRT,2))]);
    text(2-.25, walkRT*1.05, [num2str(round(walkRT,2))]);
    ylim([0 1])

    %% SDT
    subplot(2,3,5); cla;
       bar([dp_S, dp_N; crit_S, crit_N]);
    hold on;  
    title([subjID ]);
    set(gca,'xticklabels', {'dPrime', 'criterion'})
    ylabel('Rate');
    xlabel('DV')
%     ylim([0 1]);
    legend('slow', 'normal.');
    hold on;
    ylim([0 2.5])




    GFX(ippant).slow_RT = walkSRT;
    GFX(ippant).norm_RT= walkRT;

    GFX(ippant).slow_dprime = dp_S;
    GFX(ippant).norm_dprime = dp_N;

    GFX(ippant).slow_crit = crit_S;
    GFX(ippant).norm_crit = crit_N;

    %%  plot PsychFits for walking speeds:
%     % now we can extract the types correct for each.
%     %whole exp first:
%%     exprows=standingrows;
% figure()
%%
leg=[];
for id=1:2
    switch id
        case 1
            userows= walkSrows;
            col= [.2 .2 .7];
        case 2
            userows = walkNrows;
            col = [.2 .7 .2];
        
    end
%%  new version using binFitPlot
%data = [intensity';response];
stimIntensity = T.toneAmp(userows)+1; % remove 0
% stimIntensity(isinf(stimIntensity))=0;
respCorrect = T.correctResponse(userows);
dataMat = [stimIntensity,respCorrect];
gamma = .5; %^ or 0.5
figureFlag=0; % we'll plot out own.
nBins=7; 
StimList=1:7;
[ alpha , beta, Rsq , binnedData, plotOUT] = binFitPlot( dataMat, gamma, figureFlag , nBins);
%%

!TO DO:
% 
%  subplot(2,3,[3,6]);
% title('binFit (posIdx) PF');
% % axes
% hold on
% leg(id)=plot(plotOUT.xDim,plotOUT.pfit,'-','color',col,'linewidth',4);
% plot(plotOUT.xbins,plotOUT.ybins,['.-'],'color', col,'markersize',40);
% set(gca, 'fontsize',12);
% set(gca, 'Xtick',StimList, 'XTickLabelRotation',-20);
% % axis([min(StimLevels) max(StimLevels) .4 1]);
% xlabel('Stimulus Intensity');
% 
% ylabel('proportion correct');
% ylim([.0 1]);
% hold on
% % %% can we also use the contrast (raw)
% % stimIntensity = T.targContrast(userows);
% % % stimIntensity(isinf(stimIntensity))=0;
% % respCorrect = T.targCor(userows);
% % dataMat = [stimIntensity,respCorrect];
% % gamma = .5; %^ or 0.5
% % figureFlag=1; % we'll plot out own.
% % nBins=7; 
% % StimList=1:7;
% % [ alpha , beta, Rsq , binnedData, plotOUT] = binFitPlot( dataMat, gamma, figureFlag , nBins);
% % %%


end
legend(leg,{'slow', 'normal'}, 'Location','SouthEast','autoupdate','off')
%%
    print('-dpng', [subjID ' quick summary (walk speeds) '])
    
end % participant



% save('GFX_Calibration', 'GFX');
%% quick plot, similar to the bar charts above. showing the overall effects (if any). 

%panel 1: Accuracy
%panel 2: HR'
%panel 3: FAR
%panel 3: RT
%panel 4: Dprime
%panel 5: criterion

dvfields={'acc', 'HR', 'FAr', 'RT','dprime','crit'};
ymax = [1,1,1,2,3,3];

clf

fsize=15;
groupname= {'wEye', 'woEye', 'combined'};

plotGFX=GFX;
for iDV= 1:6

%wrangle the data.

tmpdata = [plotGFX.(['slow_' dvfields{iDV}]); plotGFX.(['norm_' dvfields{iDV}])]';
nppants = size(tmpdata,1);



subplot(2,3,iDV);
mData= nanmean(tmpdata,1);
errData = nanstd(tmpdata,0,1)./sqrt(nppants); % sE
bar(mData);
hold on;
errorbar(1:2, mData, errData, 'k', 'linewidth',2,'LineStyle',' none');

%ttest?
[h,p,stat] = ttest(tmpdata(:,1), tmpdata(:,2));

ylim([0 ymax(iDV)]);
ylabel(dvfields{iDV});
title(dvfields{iDV});
text(1.5, ymax(iDV)*.95, ['p = ' num2str(round(p,3))], 'HorizontalAlignment', 'center');

set(gca,'fontsize', fsize,'XtickLabels', {'slow', 'normal'});

shg
end % iDV
%
subplot(2,3,1);
title(['N(' num2str(nppants) ') ' groupname{igroup}]);
print('-dpng', ['Group Calibration summary N(' num2str(nppants) ') '  groupname{igroup} ]);

%%