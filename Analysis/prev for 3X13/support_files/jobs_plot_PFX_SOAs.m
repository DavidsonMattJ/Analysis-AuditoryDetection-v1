% jobs_plot_PFX_SOAs

% -- called from JOBS_import_preprocess;


%% THIS SECTION UPDATES FILE PATH FOR DESKTOP OR ONLINE VERSION
% 
% MatlabOnline = 0;
% 
% if MatlabOnline == 0
%     rootPath = 'C:/Users/mobil/Gabriel/MATLAB_Drive';
% else
%     rootPath = 'MATLAB DRIVE';
% end

%% PRX SOA Script begins here

% this script loads the presaved participant level data (SOA variables).
% makes a plot per participant, comparing slow and fast walk conditions.

% - load data
% - extract trials per SOA type
% - calculate proportion same
% save
%%
% load data, do x y z
% GF+

fntsize= 12;
cd(savedatadir)
allppants = dir([pwd filesep 'p_*']);
%%
for ippant= 1:length(allppants)

cd(savedatadir);
load(allppants(ippant).name, 'ppantData_SOAs', 'subjID');

% set up a figure
figure(1); 
clf;
set(gcf, 'color', 'w', 'units', 'normalized', 'position', [0 0 .6 .6]);

SOAs= ppantData_SOAs.SOAs;
nSOAs = length(ppantData_SOAs.SOAs);
% nGaitPhases = 4;
%subjID = string(subjID);

% first plot, grand average proportion same.
figure(1); clf

    hold on

plottybrah = plot(SOAs, ppantData_SOAs.propSame_all,'g-o', 'LineWidth',2);
title([ subjID ' - Same Responses Across All Walking Speeds']);
xlabel("SOA (auditory leading <-> visual leading )")
ylabel("Proportion 'Same' Response")

for i=1:nSOAs

    xlocation = SOAs(i);
    ylocation = ppantData_SOAs.propSame_all(i);
    ourText = ['n = ' num2str(ppantData_SOAs.nTrialsSOA_all(i))];
    tx = text(xlocation,ylocation,ourText);
    tx.FontWeight = "bold";
end

set(gca,'fontsize', fntsize)
ylim([0 1])
cd(figdir);
cd('SynchFunction_allPropSame')
print('-dpng', [subjID '_' num2str(ippant) '_allPropsame']);
%% start various plots
% plot slow vs fast (proportions)
figure(1)
clf
hold on
plot(SOAs,ppantData_SOAs.propSame_slow,'bo-', 'linew', 2)
plot(SOAs,ppantData_SOAs.propSame_fast,'ro-','linew', 2)
xlabel('Stimulus Onset Asynchrony (Seconds)')
ylabel('Proportion of "Same" Responses')
title([ subjID '  Same Responses both All Walking Speeds']);
legend('slow', 'fast')
set(gca,'fontsize', fntsize)
ylim([0 1])
%%

cd(figdir);
cd('SynchFunction_perSpeed')
print('-dpng', [subjID '_' num2str(ippant) '_comparePropsame']);
%% plot Response Times

figure(1)
clf
hold on
title(string(subjID) + ' - Response Time')
xlabel('Stimulus Onset Asynchrony (Seconds)')
ylabel('Response Time (Seconds)')

plot(SOAs,ppantData_SOAs.rt1_slow, 'bo-');
plot(SOAs,ppantData_SOAs.rt1_fast, 'ro-');

%plot(SOAs,meanRTSlowThisSOA, 'bo-')
%plot(SOAs,meanRTFastThisSOA, 'ro-')


legend('Slow Walk Speed', 'Fast Walk Speed')

set(gca,'fontsize', fntsize)

cd(figdir);
cd('RTs_perSOA');
print('-dpng', [subjID '_' num2str(ippant) '_ResponseTimes']);
%% Response Times adjusted for SOA - As longer SOA's require a longer time before both stimuli have been shown, this may
% lead to longer Response times as participants may wait to see both stimuli before responding.  This is accounted for by
% subtracting the SOA from the response time.

figure(1)
clf
hold on
title(string(subjID) + ' - Response Time adjusted for SOA')
xlabel('Stimulus Onset Asynchrony (Seconds)')
ylabel('Response Time minus SOA (Seconds)')
plot(SOAs,ppantData_SOAs.rt2_slow, 'bo-')
plot(SOAs,ppantData_SOAs.rt2_fast, 'ro-')
legend('Slow Walk Speed', 'Fast Walk Speed')

set(gca,'fontsize', fntsize)

cd(figdir);
cd('RTs_fromstim2');
print('-dpng', [subjID '_' num2str(ippant) 'ResponseTimesAdjustdforSOA']);

%% Prop Same by Gait Cycle and SOA

speedConds = {'_all','_slow','_fast'};
speedtitle = {'All', 'Slow', 'Normal'};
% change colormap
try cmap= brewermap(14, 'YlOrRd');
catch cmap = cbrewer('seq', 'YlOrRd', 14);
end
    cmap=cmap(6:14,:);
    %%
for icond = 2%:length(speedConds)

%% Plot PropSame by Gait Phase for each SOA

figure(1)
clf
hold on

% colormap(cmap(7:14,:));

title({[subjID ' Prop Same Responses by Gait Q & SOA'];['(Speed: ' speedtitle{icond} ')']});


dataIN = ppantData_SOAs.(['propSameByGait' speedConds{icond} '_FirstAligned']);

%%
figh=[];
for k = 1:nSOAs
    %plot(propsamebyGaitandSOA(:), propsamebyGaitandSOA(:,k), styles(k));
    %plot(propsameThisCond(:,k), styles(k));
   figh(k)= plot(1:size(dataIN,1), dataIN(:,k),'o-', 'color', cmap(k,:), 'linew', 2);

   text(0.4, dataIN(1,k), ['SOA: ' num2str(SOAs(k))]);
end
% xlim([.5 4.5]);
set(gca, 'xtick', 1:size(dataIN,1));%, 'XTickLabel', {'1', '2', '3', '4'});
xlabel('Gait phases')
ylim([0 1])
xlim([ 0 size(dataIN,1)+1])
ylabel('Proportion "Same"');

set(gca,'fontsize', fntsize)
shg
%%
legend(figh, {num2str(SOAs)});
shg
%%
%     legendinfo = SOAs;
%     legend(legendinfo);
xlabel('Gait Phase');
ylabel('Prop. "Same" responses');
hold off

filenamePropSameGatePhaseSOA = [subjID '_'  num2str(ippant)  '_PropSame_by_GaitPhase_perSOA_' speedtitle{icond}];

cd(figdir);
cd('PropSamebyGait_perSOA')
print('-dpng', filenamePropSameGatePhaseSOA);  

%% Plot PropSame by SOA for each Gait Phase

figure(1)

clf
hold on

titletext2 =[subjID  '-Prop Same Responses by SOA per Gait Q- Walk Speed:'  speedtitle{icond}];

title(titletext2);
styles = {'bo-', 'ro-', 'yo-', 'ko-', 'mo-', 'co-', 'go-'};

for iGait = 1:size(dataIN,1)
    %plot(propsamebyGaitandSOA(:), propsamebyGaitandSOA(:,k), styles(k));
    %plot(propsameThisCond(:,k), styles(k));
    plot(SOAs,dataIN(iGait,:),'-o','color', cmap(iGait,:), 'LineWidth',2);
    

end
  ylim([0 1])
legend;
%%
%     legendinfo = SOAs;
%     legend(legendinfo);
xlabel('SOA');
ylabel('Prop. "Same" responses');
hold off
set(gca,'fontsize', fntsize)
%%
cd(figdir);
cd('SynchFunction_perGaitQ')
filenamePropSameSOAGatePhase = [subjID  '_'  num2str(ippant)  '_PropSame_by_SOA_perGaitPhase-' speedtitle{icond} ];
print('-dpng', filenamePropSameGatePhaseSOA); 


end % End of Walk Speed Condition Loop

end% End of Participant Loop